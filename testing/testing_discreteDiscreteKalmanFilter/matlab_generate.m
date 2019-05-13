%script to test the equivalence between the matlab and c++ 
%implementation of the discrete discrete Extended Kalman filter

addpath('../../matlab_implementation/discrete_discrete/')

A = [-0.0004,0;0,-0.002];
aa1 = 0.01;
aa2 = 0.02;
A2 = [aa1,0;0,aa2];
next_func = @(x,t) x + A*x + 0.*t - [aa1*x(1)^2;aa2*x(2)^2];
jacobian_func = @(x,t) eye(2) + A + 0.*t + 0.*x - 2.*[aa1*x(1),0;0,aa2*x(2)]; 

R_d = 0.1^2;
Q_d = 0.01^2.*eye(2); 
C = [1,1];


state_count = 2;
sensor_count = 1;

x_0 = [1;1];
outputs = 1000;
ideal_data = zeros(state_count,outputs);
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
	ideal_data(:,ii) = x;
	process_noise_data(:,ii) = x_noise;
	x = next_func(x,0);%t=0 because function is time independent
	x_noise = next_func(x_noise,0) + sqrt(Q_d)*randn(state_count,1);
	measurements(ii) = C*x_noise + sqrt(R_d)*randn(sensor_count,1);
end

P_0 = [1,0;0,1];
t_start = 0;

[estimates, covariances] = ddekf(next_func,jacobian_func,1,t_start,state_count,...
	sensor_count,outputs,C,chol(Q_d)',chol(R_d)',chol(P_0)',x_0, measurements);

%write files for c++ code to read and run
csvwrite('Q_d.txt',Q_d);
csvwrite('R_d.txt',R_d);
csvwrite('C.txt',C);
csvwrite("qc_nonlin.txt", -[aa1,0,0;0,0,aa2] );
csvwrite("lc_nonlin.txt", eye(state_count) + A);
csvwrite("cc_nonlin.txt", zeros(state_count,1));
csvwrite("lc_jac.txt", -2*[aa1,0,0,0;0,0,0,aa2]);
csvwrite("cc_jac.txt", eye(2) + A);
csvwrite('x_0.txt',x_0);
csvwrite('P_0.txt',P_0);
csvwrite('measurements.txt',measurements);
csvwrite('measurement_count.txt',length(measurements));
csvwrite('state_count.txt',state_count);
