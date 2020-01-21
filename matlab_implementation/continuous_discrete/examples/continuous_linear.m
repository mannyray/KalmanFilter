clear all;
addpath('..');

C = [1, 0.1]; %observation matrix to obtain measurement y (y = C*x + noise)
A = [-1,0.2;-0.1,-1]; %(d/dt)x(t) = A*x(t)
R_c = [0.1^2]; %continuous sensor noise
Q_c = [0.001^2,0;0,0.001^2]; %continuous process noise
t_start = 0;
t_final = 20;
outputs = 10000;
P_0 = eye(2);%initial covariance

%same parameters as for steady_state and continuous_continuous except sensor
%matrix is scaled based on time distance between measurements 
dt = (t_final-t_start)/outputs;
R_d = (1/dt).*R_c;

%generate data and measurements 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x_0 = [100; 80]; %initial real data state
state_count = length(x_0);
sensor_count = 1;
ideal_data = zeros(state_count,outputs); %data with no process noise
process_noise_data = zeros(state_count,outputs);
measurements = zeros(sensor_count,outputs);
x = x_0;
x_noise = x_0;

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave==true
	randn("seed",1);
else
	rng(1);
end

for ii=1:outputs
	x = x + dt.*(A*x);%Explicit Eulers
	ideal_data(:,ii) = x;%the first column of ideal_data and process_noise_data is x_1
	x_noise = x_noise + dt.*(A*x_noise) + sqrt(dt.*Q_c)*randn(state_count,1);%Euler-Maruyama
	process_noise_data(:,ii) = x_noise;
	measurements(ii) = C*x_noise + sqrt(R_d)*randn(sensor_count,1);
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%set up arguments
if isOctave==true
	f_func = @(x,t) A*x + 0.*t;
	jacobian_func = @(x,t) A + 0.*x + 0.*t;
else
	f_func = @(x,t) A*x + 0.*t;
	xx = sym('x',[state_count,1]);
	jacobian_func = matlabFunction(jacobian(A*xx,xx),'Vars',{xx});
	jacobian_func = @(x,t) jacobian_func(x) + 0.*t;%no dependence on t
end

rk4_steps = 4;%type 'help cdekf' for more info 
[estimates,covariances] = cdekf(f_func,jacobian_func,dt,rk4_steps,t_start,...
	state_count,sensor_count,outputs,C,chol(Q_c)',chol(R_d)',chol(P_0)',x_0,measurements);
	
covariances{end}
