clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;
close all;

C_ukf = @(x,i) C*x + 0.*i;
Q_d_ukf = @(i) Q_d + 0.*i;
R_d_ukf = @(i) R_d + 0.*i;
func = @(x) func(x,0);
[estimates_ukf, covariances_ukf] = ukf(func, state_count, sensor_count,...
	outputs,C_ukf,Q_d_ukf,R_d_ukf,P_0,x_0, measurements);%no need for jacobian 
covariances_ukf{end}
