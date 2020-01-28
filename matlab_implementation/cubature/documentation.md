# Cubature Kalman filter

Here we will describe the Cubature Kalman filter. The algorithm is similar to the UKF, but different. The algorithm is described here

Arasaratnam, Ienkaran, and Simon Haykin. "Cubature kalman filters." IEEE Transactions on automatic control 54.6 (2009): 1254-1269.

The state evolves discretely and the measurements are discrete in time. The code is in `matlab_implementation/cubature` of [https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter). The function header in `cubature.m` is:


```Matlab
function [estimates, covariances ] = cubature(f_func,dt_between_measurements,start_time,state_count,sensor_count,measurement_count,C_func,Q_root,R_root,P_0_root,x_0, measurements)
```

For those familiar with the Kalman filter and notation are familiar with the naming of the variables. However, to be extra sure it is always best to run `help cubature`. We will break down an example below.

## Example 1

We will use the same model as in the `DD-EKF` tab - discrete logistic growth model with discrete measurements. The example is located in `matlab_implementation/unscented/examples/logistic2.m` and first runs the code from `matlab_implementation/discrete_discrete/examples/logistic.m` to define system parameters and then runs the cubature. The cubature only takes models that are state dependent and not 'time' dependent - this is reflect in the example code. The modification to make the code 'time' dependent is simple and easy way for contributing to the project.

```
clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

logistic;
close all;


C_func = @(x) C*x;

[estimates_cubature, covariances_cubature] = cubature(next_func,dt,t_start,state_count,sensor_count,...
  outputs,C_func,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
  
max_rel_err = max(abs(estimates - estimates_cubature)./estimates)
```

The error computed is `3.43968286702904e-06`.


## Example 2

We will run `Example 2` from `DD-EKF` tab in `matlab_implementation/discrete_discrete/examples/linear.m` for the `cubature` in `matlab_implementation/cubature/examples/linear2.m`:

```
clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;
close all;

[estimates_cubature, covariances_cubature] = cubature(func,dt,t_start,state_count,sensor_count,...
	outputs,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
```

to produce

```
ans =

   5.05404333754918e-07   2.45085320157428e-08
   2.45085320157428e-08   4.98056883303343e-07
```

## Example 3

In this example we use the target tracking example from the cubature paper implemented in `matlab_implementation/cubature/examples/target_tracking.m`

Arasaratnam, Ienkaran, and Simon Haykin. "Cubature kalman filters." IEEE Transactions on automatic control 54.6 (2009): 1254-1269.

and use it to compare the `UKF` and `cubature`. In the original paper, the example was only done for the square root variant of `cubature` - here it we do it for both filters with non square root variants.

We will in the figure below the estimate of the `cubature` diverges:

<img src="assets/cubature/examples/diverge1.png" alt="drawing2" style="width:80%;max-width:600px;"/>



