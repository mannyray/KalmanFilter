clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

linear;
close all;

particle_count = 30000; 
particle = mvnrnd(x_0,P_0,particle_count)'; 
[estimates_particle, covariances_particle] = particle_filter(func,jacobian_func,dt,t_start,state_count,sensor_count,...
	outputs,particle_count,particle,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
covariances_particle{end}
