# Unscented Kalman Filter

Here we will describe the Unscented Kalman filter. The state evolves discretely and the measurements are discrete in time. The code is in `matlab_implementation/unscented` of [https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter). The function header in `ukf.m` is:

```
function [estimates, covariances] = ukf(f_func, state_count, sensor_count, measurement_count,C,Q,R,P_0,x_0, measurements)
```

For those familiar with the Kalman filter and notation are familiar with the naming of the variables. However, to be extra sure it is always best to run `help ukf`. We will break down an example below.

## Example 1

We will use the same model as in the `DD-EKF` tab - discrete logistic growth model with discrete measurements. The example is located in `matlab_implementation/unscented/examples/logistic2.m` and first runs the code from `matlab_implementation/discrete_discrete/examples/logistic.m` to define system parameters and then runs the ukf. The ukf only takes models that are state dependent and not 'time' dependent - this is reflect in the example code. The modification to make the code 'time' dependent is simple and easy way for contributing to the project.

```
clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

logistic;
%close plots generated as they are not needed here.
close all;

C_ukf = @(x,i) C*x + 0.*i;
Q_d_ukf = @(i) Q_d + 0.*i;
R_d_ukf = @(i) R_d + 0.*i;
%ukf implementation here does not allow for a time argument 
%this can be implemented
next_func = @(x) next_func(x,0);

[estimates_ukf, covariances_ukf] = ukf(next_func, state_count, sensor_count,...
	 outputs,C_ukf,Q_d_ukf,R_d_ukf,P_0,x_0, measurements);
	 
max_rel_err = max(abs(estimates - estimates_ukf)./estimates)
```

The max relative error is between estimates of `ddekf` and `ukf` is `0.011142` or around one percent. `ukf2` is an equivalent implementation to `ukf` but gets rid of the clunky matrix wrappers (for example `Q_d_ukf`). The wrappers are useful if you have 'time' dependent but if your matrices are fixed and you are looking for a faster solution then go for `ukf2`. For a lot of the `Matlab` code, the code and examples are such that the noise covariance matrices and measurement matrices are usually fixed - by observing the different between `ukf2` and `ukf` you can appropriately modify the code to fit your exact purposes.

```
[estimates_ukf2, covariances_ukf2] = ukf2(next_func, state_count, sensor_count, outputs,C,Q_d,R_d,P_0,x_0, measurements);

max_err = max(abs(estimates_ukf - estimates_ukf2))
```

The error computed is `1.4211e-14`.


## Example 2

We will run `Example 2` from `DD-EKF` tab in `matlab_implementation/discrete_discrete/examples/linear.m` for the `UKF` in `matlab_implementation/unscented/examples/linear2.m`:

```
clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;

C_ukf = @(x,i) C*x + 0.*i;
Q_d_ukf = @(i) Q_d + 0.*i;
R_d_ukf = @(i) R_d + 0.*i;
func = @(x) func(x,0);
[estimates_ukf, covariances_ukf] = ukf(func, state_count, sensor_count,...
	outputs,C_ukf,Q_d_ukf,R_d_ukf,P_0,x_0, measurements);%no need for jacobian 
covariances_ukf{end}
```

to produce

```
ans =

   5.0540e-07   2.4509e-08
   2.4509e-08   4.9806e-07

ans =

   5.0540e-07   2.4509e-08
   2.4509e-08   4.9806e-07
```

The first is output by `linear` script and the second is the one computed by `ukf`.
