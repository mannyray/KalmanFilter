clear all;
addpath('..');

%setting the seed for reproducability
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave==true
	randn("seed",1);
else
	rng(1);
end

%define continuous logistic growth model and its jacobian
rate = 0.5; 
max_pop = 100;
deriv_func = @(x) rate*x*(1 - x/max_pop);%logistic population model (quadratic in x)
if isOctave==true
	jacobian_func = @(x) rate - (2*x*rate)/max_pop;%manual derivative - octave might have its own function
else
	%in case you don't want to compute the jacobian manually and you have access
	%to Matlab then you can compute it symbolically and then extract a function from it
	xx = sym('x',[1,1]);
	jacobian_func = matlabFunction(jacobian(deriv_func(xx),xx),'Vars',xx);
end
deriv_func = @(x,t) deriv_func(x) + 0.*t; 
jacobian_func = @(x,t) jacobian_func(x) + 0.*t;
 

%in this model the time is only necessary
%for plotting an interpreting the results.
t_start = 0;
t_final = 10;
outputs = 1000;
dt = (t_final - t_start)/outputs;
 

R_c = 0.01; 
%discrete noise covariance. 
R_d = 1;
%continuous process noise covariance
Q_c = 10;
%measurement matrix
C = 1;
 
%generate data and measurements
state_count = 1; 
sensor_count = 1;

%model results if there was no process noise
ideal_data = zeros(state_count,outputs);
%actual data (with process noise) 
process_noise_data = zeros(state_count,outputs);
%noisy measurements of actual data
measurements = zeros(sensor_count,outputs);

%initial condition for ideal
%and real process data
x_0 = max_pop/2;
P_0 = 1;

%initial condition for ideal
%and real process data
x = x_0;
x_noise = x_0;

%generate the data
curTime = t_start;
for ii=1:outputs
	%Explicit Eulers
	x = x + dt.*deriv_func(x,curTime);
	ideal_data(:,ii) = x;
	
	%Euler-Maruyama
	x_noise = x_noise+dt.*deriv_func(x_noise,curTime)+sqrt(dt.*Q_c)*randn(state_count,1);
	process_noise_data(:,ii) = x_noise;
	measurements(ii) = C*x_noise + sqrt(R_d)*randn(sensor_count,1);
	
	curTime = curTime + dt;
end

%filter the noisy measurements
rk4_steps = 1;%internal for code - type 'help cdekf' for more info 
[estimates,covariances] = cdekf(deriv_func,jacobian_func,dt,rk4_steps,t_start,...
	state_count,sensor_count,outputs,C,chol(Q_c)',chol(R_d)',chol(P_0)',x_0,measurements);

%The ideal data is not plotted and is only for reference
time = 0:dt:(outputs-1)*dt;
h = figure;
hold on;
limit = 100;
plot(time(1:limit),measurements(1:limit),'LineWidth',2); 
plot(time(1:limit),process_noise_data(1:limit),'LineWidth',2);
plot(time(1:limit),estimates(1:limit),'LineWidth',2);
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')

h = figure;
hold on;
plot(0:dt:(outputs-1)*dt,measurements); 
plot(0:dt:(outputs-1)*dt,process_noise_data);
plot(0:dt:(outputs-1)*dt,estimates(1:end-1));
legend('Measurement','Real data','Estimate');
xlabel('Time')
ylabel('Population')
