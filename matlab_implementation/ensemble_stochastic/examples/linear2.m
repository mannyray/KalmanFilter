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
