clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

logistic;
close all;

C_ukf = @(x,i) C*x + 0.*i;
Q_d_ukf = @(i) Q_d + 0.*i;
R_d_ukf = @(i) R_d + 0.*i;
next_func = @(x) next_func(x,0);

[estimates_ukf, covariances_ukf] = ukf(next_func, state_count, sensor_count,...
	 outputs,C_ukf,Q_d_ukf,R_d_ukf,P_0,x_0, measurements);
	 
max_err = max(abs(estimates - estimates_ukf)./estimates)

[estimates_ukf2, covariances_ukf2] = ukf2(next_func, state_count, sensor_count, outputs,C,Q_d,R_d,P_0,x_0, measurements);

max_err = max(abs(estimates_ukf - estimates_ukf2))
