# Ensemble Stochastic

Here we will describe the stochastic Ensemble Kalman Filter. The state evolves discretely and the measurements are discrete in time. The code is in `matlab_implementation/ensemble_stochastic` of [https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter). The function header in `enkf_stochastic.m` is:

```
function [estimates,covariances] = enkf_stochastic(func,state_count,sensor_count,measurement_count,ensemble_size,C,R,Q,P_0,x_0, measurements)
```

For those familiar with the Kalman filter and notation are familiar with the naming of the variables. However, to be extra sure it is always best to run `help enfk_stochastic`. We will break down an example below. The `enkf` only takes models that are state dependent and not 'time' dependent. The modification to make it 'time' dependent is simple and easy to do for contributing to the project.
 
## Example 1

We will use the same model as in the `DD-EKF` tab. The example is located in `matlab_implementation/ensemble_stochastic/examples/logistic2.m` and first runs the code from `matlab_implementation/discrete_discrete/examples/logistic.m` to define system parameters and then runs the ensemble filter. Run `matlab_implementation/ensemble_stochastic/examples/logistic2.m`:

```
logistic;
close all;

ensemble_size = 5;
next_func = @(x) next_func(x,0);
[estimates_enkf,covariances_enkf] = enkf_stochastic(next_func,state_count,sensor_count,...
	outputs,ensemble_size,C,R_d,Q_d,P_0,x_0, measurements);
```

Now plot the results

```
h = figure;
hold on;
plot(0:dt:(100-1)*dt,measurements(1:100),'LineWidth',2); 
plot(0:dt:(100-1)*dt,process_noise_data(1:100),'LineWidth',2);
plot(0:dt:(100-1)*dt,estimates_enkf(1:100),'LineWidth',2);
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')
```
<img src="../../ensemble_stochastic/examples/logistic1.png" alt="drawing2" style="width:80%;max-width:600px;"/>






## Example 2

We will run `Example 2` from `DD-EKF` tab in `matlab_implementation/discrete_discrete/examples/linear.m` for the `enkf_stochastic` filter in `matlab_implementation/ensemble_stochastic/examples/linear2.m`:

```
clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;
close all;

func = @(x) func(x,0);
ensemble_size = 100000;
[estimates_enkf,covariances_enkf] = enkf_stochastic(func,state_count,sensor_count,...
	outputs,ensemble_size,C,R_d,Q_d,P_0,x_0, measurements);
covariances_enkf{end}
```

to produce

```
ans =

   5.0540e-07   2.4509e-08
   2.4509e-08   4.9806e-07
```

The first is output by `linear` script and the second is the one computed by ensemble filter.
