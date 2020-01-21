clear all;
addpath('..');
addpath('../../discrete_discrete');
addpath('../../discrete_discrete/examples');

logistic;
close all;

ensemble_size = 5;
next_func = @(x) next_func(x,0);

[estimates_enkf,covariances_enkf] = enkf_stochastic(next_func,state_count,sensor_count,...
	outputs,ensemble_size,C,R_d,Q_d,P_0,x_0, measurements);

h = figure;
hold on;
plot(0:dt:(100-1)*dt,measurements(1:100),'LineWidth',2); 
plot(0:dt:(100-1)*dt,process_noise_data(1:100),'LineWidth',2);
plot(0:dt:(100-1)*dt,estimates_enkf(1:100),'LineWidth',2);
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')
