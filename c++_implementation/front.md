# Kalman Filter

This is a `C++` Kalman filter library designed to work with various data types and provide flexibility to use more advanced models. The `C++` language may seem verbose to those transitioning from say `Python` or `Matlab`. However, this library may seem verbose not only to the `Python` user but to the `C++` user as well. To address this concern, many examples, good documentation and motivation behind the code will be discussed. This is an active project and will continue to grow. I am planning to add additional filters such as the UKF and add further examples. I am looking for help in optimizing the code further as well as adding additional filters and features - if you have ideas, advice or bugs then let me know.


If you are new to Kalman filters then please check out rlabbe's [github repository](https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python) to get an intuition behind the Kalman filter. I assume if you are reading past this point, you have a basic understanding of the Kalman filter, what it requires and common notation used.


First thing let's install the code. The code has been tested on Ubuntu system for now. 

```
$ git clone https://github.com/mannyray/KalmanFilter.git
$ cd KalmanFilter/c++_implementation
$ chmod u+x install.sh
$ ./install.sh
```

The install script only modifies `include/Eigen.h`. This library is a header-only library.

We will be looking at `examples/basicExample.cpp` and `examples/sampleModel.h` in the following tutorial to understand the library. The Kalman filter can be described by the following image (page 147) from `Introduction to Random Signals and Applied Kalman Filtering with Matlab Exercises` by `Robert Grover Brown` and `Patrick Y. C. Hwang`:


![](filter1.png)

To run this we need a few ingredients. Using the notation from the image, first we need the _transition equation_ for the 'Project ahead' block:

\f$ \hat{x}_{k+1}^- = \phi_k \hat{x}_{k} \f$  .

In addition, from the 'Project ahead' we need the _transition Jacobian_ \f$ \bf{\Phi}_k \f$ and _process noise_ \f$ \bf{Q}_k \f$. We need a _measurement equation_ that is not included in the diagram here:

\f$ z_k = C(x_k) \f$   .

Finally, we need the _measurement Jacobian_ \f$ H_k \f$ and sensor noise \f$R_k\f$. Let's consider a specific system:

\f[
	\ddot{x} + \sqrt{2}w_0\dot{x} + w_0^2x = b\sigma u
\f]

This ODE describes a damped harmonic oscillator at position \f$ x \f$ with a forcing white noise term on the right hand side. The equation is taken from `Introduction to Random Signals and Applied Kalman Filtering with Matlab Exercises` section 4.4 `Robert Grover Brown` and `Patrick Y. C. Hwang` book mentioned earlier. The model can be expressed as follows
\f[
	\frac{d}{dt}\begin{bmatrix}x_1 \\ x_2\end{bmatrix} = \begin{bmatrix} 0 & 1\\ -w_0^2 & -\sqrt{2}w_0 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
\f]
where \f$x_1 = x \f$ and \f$ x_2 = \dot{x} \f$. Using the Van Loan method described in section 3.9 of the afformentioned source the model can be discretized assuming \f$ \Delta t = 1\f$ using the following Octave script:

~~~~~~~~~~~~~~~{.m}
pkg load symbolic
syms w_0 
syms b % b = sqrt(2 * sqrt(2) w_0^3)
syms sig
syms delta_t
F = [sym(0),1;-w_0*w_0, -sym(14142135623731)/10000000000000*w_0];
G = [0;b*sig];
A = delta_t*[-F,G*G';[0,0;0,0],F'];
B = expm(A);
delta_t_val = sym(1)/1;
w_0_val = sym(1)/10;%use whatever values you want 
b_val = sym(1000000000000)/18803015465432;%octave requires such measures to reduces numerical errors
sig_val = sym(1)/1;
B_sub = eval(subs(subs(subs(subs(B,b,b_val),w_0,w_0_val),sig,sig_val),delta_t,delta_t_val));
transitionMatrix = B_sub(3:4,3:4)';
processNoise = real(transitionMatrix*B_sub(1:2,3:4));
~~~~~~~~~~~~~~~

to produce

\f[
	\begin{bmatrix} X_1 \\ X_2 \end{bmatrix}_{k+1} =   \begin{bmatrix} 0.9952315 & 0.9309552 \\ -0.0093096 & 0.8635746 \end{bmatrix} \begin{bmatrix}X_1 \\ X_2 \end{bmatrix}_k
\f]

where \f$(X_1)_k \approx x_1(t = k) \f$ and \f$(X_2)_k \approx x_2(t = k) \f$ with process noise

\f[
	\bf{Q}_k = \begin{bmatrix}  8.4741e-04 &  1.2257e-03 \\    1.2257e-03  & 2.4557e-03 \end{bmatrix}.
\f]
We will select to measure \f$X_1\f$ only so our measurement equation will be 
\f[
	z_k = \begin{bmatrix}1 & 0 \end{bmatrix} \begin{bmatrix}X_1 \\ X_2 \end{bmatrix}.
\f]
The sensor noise will be 
\f[
	R_k = (0.5)^2 = 0.25.
\f]
The transition jacobian will naturally be
\f[
\begin{bmatrix} 0.9952315 & 0.9309552 \\ -0.0093096 & 0.8635746 \end{bmatrix}
\f]
with measurement jacobian
\f[
\begin{bmatrix}1 & 0 \end{bmatrix}.
\f]


By now you are familiar with what the Kalman filter needs:
 - transition equation
 - measurement equation
 - process noise covariance
 - sensor noise covariance 
 - transition jacobian 
 - measurement jacobian 
 
Let's implement the ingredients for our specific example. In this case we will use the [`Eigen`](http://eigen.tuxfamily.org/index.php?title=Main_Page) library. 
 
### Transition equation

All ingredients or components of a system will inherit from general ingredient or component classes. In this example, we will call the class `stateModel` and it will inherit from template class `discreteModel` of `Eigen::VectorXd` type since we are working with the `Eigen` library. The general class `discreteModel` is presented first:

~~~~~~~~~~~~~~~{.cpp}
template <typename T> class discreteModel: public model<int>{
        public:
                virtual T function(const T& val, const int index) const = 0;
};
~~~~~~~~~~~~~~~
`discreteModel` is abstract and any class that inherits from it needs to implement the `function` method. We can use this class to easily implement (1):

~~~~~~~~~~~~~~~{.cpp}
class stateModel: public discreteModel<Eigen::VectorXd>{
	Eigen::MatrixXd v;
	public:
		stateModel(){
			Eigen::MatrixXd tmp(2,2);
			tmp<<0.995231,0.9309552,-0.0093096,0.8635746;
			v = tmp;
		}
		Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
			return v*val;
		}
};
~~~~~~~~~~~~~~~

`stateModel` implements `function` of `discreteModel` but now has an extra member and even a constructor. This is where the flexibility to the user comes into play. If you want to pass parameters to your constructor in order to have different and multiple instances of a model then you are free to do so. Why have multiple instances? Maybe you are working on tuning your filter or are experimenting with assuming various models. No matter the reason, this library gives you the opportunity to do so easily.

`function` in this case implements the important relationship of (1). The first argument of `function` within the Kalman filter will be the previous estimate of the model while `time` is the current model time which allows for an implementation of a time based model.

### Measurement equation

In this instance we implement (2). We will call the class `measurementModel`: 

~~~~~~~~~~~~~~~{.cpp}
class measurementModel: public discreteModel<Eigen::VectorXd>{
	Eigen::MatrixXd v;
	public:
		measurementModel(){
			Eigen::MatrixXd tmp(1,2);tmp<<1.0,0.0;
			v = tmp;
		}
		Eigen::VectorXd function(const Eigen::VectorXd & val, const int time) const override{
			return v*val;
		}
};
~~~~~~~~~~~~~~~

`measurementModel` inherits from the same `discreteModel` class as did `stateModel` of the transition equation.

### Process Noise

`discreteModel` was used for (1) and (2). Similarly for both sensor and process noise the following base class will be used:

~~~~~~~~~~~~~~~{.cpp}
template <class VECTOR, class MATRIX>
class discreteNoiseCovariance: public noiseCovariance<VECTOR,MATRIX,int>{
	public:
		virtual MATRIX function(const VECTOR & val, int t) = 0; 
		virtual MATRIX sqrt(const VECTOR & val, int t) = 0;
};
~~~~~~~~~~~~~~~

The two pure virtual methods are `function` and `sqrt`. `function` returns the covariance matrix at time `t` while - the `VECTOR` input is available for use in case the covariance depends on the state/estimate of the system - such as in distance based applications. The second method is `sqrt` which will return the Cholesky decomposition  (specify which one - upper or lower) which is used for sampling noises in `discreteDiscreteFilterSolver` (see [here](https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Drawing_values_from_the_distribution) for details) and is used in the square root filter.

The process noise (3) is implemented in class `processNoise`:

~~~~~~~~~~~~~~~{.cpp}
class processNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd v;
	Eigen::MatrixXd sqrtMat;
	public:
		processNoise(){
			Eigen::MatrixXd tmp(2,2);
			tmp<<8.4741e-04,1.2257e-03,1.2257e-03,2.4557e-03; 
			v = tmp;
			Eigen::LLT<Eigen::MatrixXd> lltOfV(v);
			sqrtMat = lltOfV.matrixL();
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
			return v;
		}
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return  sqrtMat; 
		}
};
~~~~~~~~~~~~~~~

In these cases the methods are simple and return an already precomputed result without making use of the arguments as the system is simple. Despite the simplicity, it is the _user's responsibility_ to make sure the dimensions of various matrices, vectors that are returned or used in the various 'ingredients' of the system are logically correct and consistent.

### Sensor Noise

~~~~~~~~~~~~~~~{.cpp}
class sensorNoise: public discreteNoiseCovariance<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd v;
	Eigen::MatrixXd sqrtMat;
	public:
		sensorNoise(){
			Eigen::MatrixXd tmp(1,1);
			tmp<<0.25;
			v = tmp;
			Eigen::MatrixXd tmp2(1,1);
			tmp2<<0.5;
			sqrtMat = tmp2; 
		}
		Eigen::MatrixXd function(const Eigen::VectorXd &est, int t) override{
			return v;
		}
		Eigen::MatrixXd sqrt(const Eigen::VectorXd &est, int t) override{
			return sqrtMat; 
		}
};
~~~~~~~~~~~~~~~

### Transition Jacobian
~~~~~~~~~~~~~~~{.cpp}
class transitionJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd v;
	public:
		transitionJac(){
			Eigen::MatrixXd tmp(2,2);
			tmp<<0.995231,0.9309552,-0.0093096,0.8635746;
			v = tmp;
		}

		Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
			return v;
		}
};
~~~~~~~~~~~~~~~

### Measurement Jacobian

~~~~~~~~~~~~~~~{.cpp}
class measurementJac: public jacobianDiscrete<Eigen::VectorXd,Eigen::MatrixXd>{
	Eigen::MatrixXd v;
	public:
		measurementJac(){
			Eigen::MatrixXd tmp(1,2);
			tmp<<1.0,0.0;
			v = tmp;
		}

		Eigen::MatrixXd function(const Eigen::VectorXd & val, int t){
			return v;
		}
};
~~~~~~~~~~~~~~~

## Putting it all together

Since we have all the ingredients (in `examples/sampleModel.h`), we are ready to combine them together to define the model. First create instances of the previously discussed classes:

~~~~~~~~~~~~~~~{.cpp}
stateModel tm;
measurementModel mm;
processNoise pn;
sensorNoise sn;
transitionJacobian tj;
measurementJacobian mj;
~~~~~~~~~~~~~~~
Then specify the state count and sensor count for this specific model

~~~~~~~~~~~~~~~{.cpp}
int stateCount = 2;
int sensorCount = 1;
~~~~~~~~~~~~~~~

Now we create an instance of the model.  

~~~~~~~~~~~~~~~{.cpp}
discreteDiscreteFilterModel<Eigen::VectorXd,Eigen::MatrixXd> ddfm(&tm,&mm,&pn,&sn,&tj,&mj,stateCount,sensorCount);
~~~~~~~~~~~~~~~

`discreteDiscreteFilterModel`'s purpose is to provide a nice collection of the ingredients in one class while providing an interface layer of the ingredient's methods to other classes. The interface layer allows for easier refactoring of the underlying ingredients.

As probably guessed by the reader, the parameters of `discreteDiscreteFilterModel`'s constructor are pointers to abstract types such as `discreteModel` which makes the code easy to use for various models/jacobians/noises. The `<Eigen::VectorXd,Eigen::MatrixXd>` specifies which vector and matrix type is used in the code as the majority of the code is built on templates - read the documentation to make sure the types and function parameters entered are in correct order.

Two classes that depend on the `discreteDiscreteFilterModel` are `discreteDiscreteFilterSolver` and the filter itself: `discreteDiscreteKalmanFilter`. `discreteDiscreteFilterSolver` is used in cases when one does not have any sample data/measurements and needs to generate them as is common when running simulations. Let's generate some data for our specific model. We need to specify an initial condition and initial time for the system as well:

~~~~~~~~~~~~~~~{.cpp}
int initialTime = 0;
Eigen::VectorXd initialState(2);
initialState<<0,0;
discreteDiscreteFilterSolver<Eigen::VectorXd,Eigen::MatrixXd> ddfs(&ddfm,initialTime,initialState);
~~~~~~~~~~~~~~~
Let's say we want to solve for `10000` time steps:
~~~~~~~~~~~~~~~{.cpp}
int steps = 10000;
ddfs.solve(steps);
~~~~~~~~~~~~~~~
To get generated measurements and the true state of the system run:
~~~~~~~~~~~~~~~{.cpp}
Eigen::VectorXd * states = ddfs.getSolvedStates();
Eigen::VectorXd * measurements = ddfs.getSolvedMeasurements()
~~~~~~~~~~~~~~~

We now have almost everything. Before running the Kalman filter we need to define our initial estimate and initial covariance:
~~~~~~~~~~~~~~~{.cpp}
Eigen::VectorXd initiaEstimate(2);
initialEstimate<<0,0;
Eigen::MatrixXd initialCovariance(2,2);
initialCovariance<<1,0,0,0.01;
~~~~~~~~~~~~~~~
In our idealized case here the `initialEstimate` and `initialState` match up exactly.

If you have read the afformentioned rlabbe repository or studied the Kalman filter previously you will know that the filter consists of the `predict` and `update` phase that are repeated for all the measurements over and over again. Hence to filter through the measurements, we run:

~~~~~~~~~~~~~~~{.cpp}
discreteDiscreteKalmanFilter<Eigen::VectorXd,Eigen::MatrixXd> filter(initialTime, initialEstimate, initialCovariance,&ddfm);

Eigen::VectorXd * estimates = new Eigen::VectorXd[steps];
for(int i = 0; i < steps; i++){
	filter.predict(1);
	filter.update(measurements[i]);
	
	estimates[i] = filter.getCurrentEstimate();
	//save the estimate or save the covariance, save only every 10th estimate,do nothing,
	//just run predict if you're simulating infrequent measurements..it's up to you and your context
}
~~~~~~~~~~~~~~~


We have finally filtered the data and the estimates are stored in `estimates` array. The user may be wondering why put the burden of running the predict and update on the user and not the Kalman Filter class. Why not simply pass in the `measurements` array into some method? First of all, if you have gotten all the way here after reviewing the ingridients and are upset with a simple for loop then shame on you. Second of all, the context in which the Kalman filter can be used is diverse. Users may have different output requirements. In addition, passing in a measurement array implies post processing of data and does not allow for 'real time' processing and logic such as adding system input. Finally, it gets the user used to the update/predict terminogy as is common with the Kalman filter and lets them avoid off-by-one logic errors.

Run `cd examples; g++ basixExample.cpp; ./a.out` to save the data into `estimates.txt`, `states.txt`, `measurements.txt`. We can plot the results using Octave:
~~~~~~~~~~~~~~~{.m}
dd = load('estimates.txt');
estimates = [dd(1:2:end)';dd(2:2:end)'];
dd = load('states.txt');
states = [dd(1:2:end)';dd(2:2:end)'];
dd = load('measurements.txt');
measurements = dd';
limits = 200;
hold on
plot(1:limits,measurements(1:limits),'r','linewidth',1);
plot(1:limits,estimates(1,1:limits),'b','linewidth',2);
plot(1:limits,states(1,1:limits),'g','linewidth',2);
xlabel ("Time");
ylabel ("Position");
legend("Measurement","Estimate","True State");
~~~~~~~~~~~~~~~
to obtain

<img src="filter2.png" alt="drawing" width="80%"/>


Let's compare the results to the matlab implementation 
~~~~~~~~~~~~~~~{.m}
addpath('../../../KalmanFilter/matlab_implementation/discrete_discrete/')
next_func = @(x,t) transitionMatrix*x;
jacobian_func = @(x,t) transitionMatrix;
dt = 1;
t_start = 0;
state_count = 2;
sensor_count = 1;
outputs = 10000;
C = [1 0];
Q_d = processNoise;
R_d = 0.25;
P_0 = [1,0;0,0.01];
x_0 = [0;0];

[estimatesOctave, covariancesOctave ] = ddekf(next_func,jacobian_func,dt,t_start,state_count,...
	sensor_count,outputs,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
~~~~~~~~~~~~~~~

Testing against the matlab square root implementation

~~~~~~~~~~~~~~~{.m}
hold on;
plot(estimates(1,1:limits),'r','linewidth',3)
plot(estimatesOctave(1,2:limits+1),'b')
plot(estimates(1,1:limits) - estimatesOctave(1,2:limits+1),'g')
xlabel ("Time");
ylabel ("Position");
legend("C++ position","Matlab square root estimate","Difference")
max(abs(estimates(1,:) - estimatesOctave(1,2:end)))%max error
max(abs(estimates(1,:) - estimatesOctave(1,2:end))./estimates(1,:))%max absolute error
sum(abs(abs(estimates(1,:) - estimatesOctave(1,2:end))./estimates(1,:)))/10000
~~~~~~~~~~~~~~~


<img src="filter3.png" alt="drawing2" width="80%"/>

with a maximum error of `1.64246433278281e-04` with max relative error of `56` percent (since they were close to zero estimates with one being `1.56686514969348e-04` while other`1.00210021201516e-04`). The average relative error is `0.11` percent. The test is passed.

### Discussing the templates

Why use templates? The idea is to allow for the use of various matrix/vector libraries in combination with a Kalman filter library. Two main ones that come to mind are `Eigen` and `Boost`. The code tutorial above is in `examples/basicExample.cpp` for the `Eigen` library while identical boost copy is located within `examples/basicExampleBoost.cpp`. After compiling and running, the results are identical:

~~~~~~~~~~~~~~~{.m}
dd = load('estimatesBoost.txt');
estimatesBoost = [dd(1:2:end)';dd(2:2:end)'];
dd = load('estimates.txt');
estimates = [dd(1:2:end)';dd(2:2:end)'];
max(estimates(1,:) - estimatesBoost(1,:))
~~~~~~~~~~~~~~~

with max difference `3.81361608958741e-14`. 

`examples/basicExamplesBoost.cpp` uses a wrapper class from `include/mathWrapper` that wraps the `boost` library. The reason for wrapping is because the template class in `include/KalmanFitler.h` uses operators for matrices and vectors such as `+` and `*` that are defined for `Eigen` but not `boost` (see [1](https://stackoverflow.com/questions/20056774/why-has-ublas-no-operatormatrix-vector), [2](http://boost.2283326.n4.nabble.com/Patch-proposal-for-overloading-operator-in-ublas-td2711622.html), [3](https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/numeric/ublas/doc/overview.htm)) or potentially other libraries. In fact, the original example `examples/basicExample.cpp` uses an extension to the `Eigen` library in `include/Eigen.h` instead of the regular `<Eigen/Dense>` since an additional random vector function needed to be added used in `discreteDiscreteFilterSolver`. The wrappers will provide the opportunity for further functionality required to implement square root versions of filters.

Instead of using `Eigen` or `Boost` lets use the good old `double` to introduce a basic one dimensional nonlinear filtering example.

### Nonlinear example:

A nonlinear example is included `examples/basicExampleNonlinear.cpp` that uses the `double` wrapper in `include/mathWrapper/double.h`. This wrapper allows for simple one dimensional examples without needing to install any extra libraries such as `Eigen` or `Boost`. The model used is the discrete logistic population growth. 

~~~~~~~~~~~~~~~{.m}
dd = load('estimatesNonlin.txt');
estimates = dd';
dd = load('statesNonlin.txt');
states = dd';
dd = load('measurementsNonlin.txt');
measurements = dd';
limits = 100;
hold on
plot(1:limits,measurements(1:limits),'r','linewidth',1);
plot(1:limits,estimates(1:limits),'b','linewidth',2);
plot(1:limits,states(1:limits),'g','linewidth',2);
xlabel ("Time");
ylabel ("Population");
legend("Measurement","Estimate","True State");
~~~~~~~~~~~~~~~
 
<img src="filter4.png" alt="drawing2" width="80%"/>
