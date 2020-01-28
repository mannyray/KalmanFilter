clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;
close all;

C_func = @(x) C*x;
[estimates_cubature, covariances_cubature] = cubature(func,dt,t_start,state_count,sensor_count,...
	outputs,C_func,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
