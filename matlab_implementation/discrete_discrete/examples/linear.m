%LINEAR EXAMPLE
clear all;
addpath('..')
C = [1, 0.1]; %observation matrix to obtain measurement y (y = C*x + noise)
A = [-1,0.2;-0.1,-1]; %(d/dt)x(t) = A*x(t)
R_c = [0.1^2]; %continuous sensor noise
Q_c = [0.001^2,0;0,0.001^2]; %continuous process noise
P_0 = eye(2);%initial covariance
state_count = 2; 

t_start = 0;
t_final = 20;
outputs = 10000;
dt = (t_final-t_start)/outputs;

A_d = eye(2) + dt.*A;%x_{k+1} = A_d*x_k
Q_d = dt.*Q_c;
R_d = (1/dt).*R_c;

func = @(x) A_d*x;
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave==true
	jacobian_func = @(x) A_d + 0.*x;
else
	xx = sym('x',[state_count,1]);
	jacobian_func = matlabFunction(jacobian(func(xx),xx),'Vars',xx);
end
func = @(x,t) func(x) + 0.*t; 
jacobian_func = @(x,t) jacobian_func(x) + 0.*t;

x_0 = [100; 80]; %initial real data state
state_count = length(x_0);
sensor_count = 1;
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
	x = x + dt.*(A*x);%Explicit Eulers
	ideal_data(:,ii) = x;%the first column of ideal_data and process_noise_data is x_1
	x_noise = x_noise + dt.*(A*x_noise) + sqrt(Q_d)*randn(state_count,1);%Euler-Maruyama
	process_noise_data(:,ii) = x_noise;
	measurements(ii) = C*x_noise + sqrt(R_d)*randn(sensor_count,1);
end

[estimates, covariances ] = ddekf(func,jacobian_func,dt,t_start,state_count,...
	sensor_count,outputs,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);
	
	
covariances{end}

