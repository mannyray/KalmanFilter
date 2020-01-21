clear all;
%include code 
addpath('..')

%setting the seed for reproducability 
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave==true
	randn("seed",1);
else
	rng(1);
end

%define discrete logistic growth model and its jacobian
rate = 0.01;
max_pop = 100;
next_func = @(x,t) x + rate.*x.*(1 - x/max_pop) + 0.*t;%time independent
jacobian_func = @(x,t) 1 + rate - (2*rate*x)/max_pop  + 0.*t;

%in this model the time is only necessary
%for plotting an interpreting the results.
outputs = 1000;
t_start = 0;
t_final = 10;
dt = (t_final - t_start)/outputs;

%discrete noise covariance
Q_d = 1;
%sensor noise covariance
R_d = 3;
%measurement matrix
C = 1;

state_count = 1; 
sensor_count = 1;

%model results if there was no process noise
ideal_data = zeros(state_count,outputs);
%actual data (with process noise) 
process_noise_data = zeros(state_count,outputs);
%noisy measurements of actual data
measurements = zeros(sensor_count,outputs);

%initial estimate and covariance
%the initial estimate in this case
%matches the actual data
x_0 = max_pop/2;
P_0 = 1;

%initial condition for ideal
%and real process data
x = x_0;
x_noise = x_0;

%generate the data
for ii=1:outputs
	x = next_func(x,ii);
	ideal_data(:,ii) = x;

	x_noise = next_func(x_noise,0) + chol(Q_d)'*randn(state_count,1);
	process_noise_data(:,ii) = x_noise;
	measurements(ii) = C*x_noise + chol(R_d)'*randn(sensor_count,1);
end

%filter the noisy measurements
[estimates, covariances] = ddekf(next_func,jacobian_func,dt,t_start,state_count,...
	sensor_count,outputs,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);

%compare the measurements, estimates and true state.
%The ideal data is not plotted and is only for reference.
time = 0:dt:(outputs-1)*dt;
h = figure;
hold on;
plot(time,measurements); 
plot(time,process_noise_data);
plot(time,estimates(1:end-1));
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')

h = figure;
hold on;
limit = 100;
plot(time(1:limit),measurements(1:limit),'LineWidth',2); 
plot(time(1:limit),process_noise_data(1:limit),'LineWidth',2);
plot(time(1:limit),estimates(1:limit),'LineWidth',2);
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')
