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
