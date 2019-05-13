clear all;
addpath('../../matlab_implementation/continuous_discrete/')
rate = 0.5; 
max_pop = 100;
deriv_func = @(x) rate*x*(1 - x/max_pop);%logistic population model (quadratic in x)

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave==true
	jacobian_func = @(x) rate - (2*x*rate)/max_pop;%manual derivative - octave might have its own function
else
	xx = sym('x',[1,1]);
	jacobian_func = matlabFunction(jacobian(deriv_func(xx),xx),'Vars',xx);
end

deriv_func = @(x,t) deriv_func(x) + 0.*t; 
jacobian_func = @(x,t) jacobian_func(x) + 0.*t;
 
t_start = 0;
t_final = 10;
outputs = 1000;
dt = (t_final - t_start)/outputs;
 
R_c = 0.01; 
R_d = (1/dt)*R_c;
Q_c = 10;
C = 1;
 
%generate data and measurements
state_count = 1; 
sensor_count = 1;
x_0 = max_pop/2;
ideal_data = zeros(state_count,outputs); %data with no process noise
process_noise_data = zeros(state_count,outputs);
measurements = zeros(sensor_count,outputs);
x = x_0;
x_noise = x_0;
 
if isOctave==true
	randn("seed",1);
else
	rng(1);
end

for ii=1:outputs
	x_noise = x_noise+dt.*deriv_func(x_noise,0)+sqrt(dt.*Q_c)*randn(state_count,1);%Euler-Maruyama
	process_noise_data(:,ii) = x_noise;
	x = x + dt.*deriv_func(x,0);%Explicit Eulers /time independent
	ideal_data(:,ii) = x;
	measurements(ii) = C*x_noise + sqrt(R_d)*randn(sensor_count,1);
end

P_0 = 1;
rk4_steps = 1;%internal for code - type 'help cdekf' for more info 
[estimates,covariances] = cdekf(deriv_func,jacobian_func,dt,rk4_steps,t_start,...
	state_count,sensor_count,outputs,C,chol(Q_c)',chol(R_d)',chol(P_0)',x_0,measurements);


csvwrite('measurements.txt',measurements);

