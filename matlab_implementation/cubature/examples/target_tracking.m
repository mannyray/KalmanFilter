clear all;
close all;
addpath('..');

%see paper 
%Arasaratnam, Ienkaran, and Simon Haykin. "Cubature kalman filters." IEEE Transactions on automatic control 54.6 (2009): 1254-1269.
%for further details

%parameters
mult_constant = pi/180;
unknown_turn_rate = -3 *mult_constant;
dt_between_measurements = 1; 
q_1 = 0.1;
q_2 = 1.75*10^(-4);
M = [dt_between_measurements^3/3, dt_between_measurements^2/2;...
	dt_between_measurements^2/2, dt_between_measurements];
	
%process noise
Q = zeros(5,5);
Q(1:2,1:2) = q_1*M;
Q(3:4,3:4) = q_2*M;
Q(5,5) = q_2*dt_between_measurements;

%sensor noise
sigma_radius = 10;
sigma_theta = sqrt(10)*0.001;
R = diag([sigma_radius^2,sigma_theta^2]);

%measurement of radius and angle
C = @(x) [sqrt(x(1)^2 + x(3)^2); atan2(x(3),x(1))];

%transition function
next_func = @(x,k)  [1,sin(mult_constant*unknown_turn_rate*dt_between_measurements)/unknown_turn_rate,0,-( (1-cos(mult_constant*unknown_turn_rate*dt_between_measurements))/unknown_turn_rate),0;...
	0, cos(mult_constant*unknown_turn_rate*dt_between_measurements), 0, -sin(mult_constant*unknown_turn_rate*dt_between_measurements),0;...
	0, (1-cos(mult_constant*unknown_turn_rate*dt_between_measurements))/unknown_turn_rate, 1, sin(mult_constant*unknown_turn_rate*dt_between_measurements)/unknown_turn_rate,0;...
	0, sin(mult_constant*unknown_turn_rate*dt_between_measurements),0,cos(mult_constant*unknown_turn_rate*dt_between_measurements),0;...
	0,0,0,0,1]*x;

%initial estimate and covariance
x_0 = [1000;300;1000;0;-3*mult_constant];
P_0 = diag([100,10,100,10,100*0.001]);

outputs = 5000;

ideal_x=x_0;
ideal_data = zeros(5,outputs);
real_x=x_0;
real_data = zeros(5,outputs);
time = 0;
measurements = zeros(2,outputs);
for i=1:outputs
	time = time + dt_between_measurements;
	ideal_x = next_func(ideal_x,time);
	ideal_data(:,i) = ideal_x;
	real_x = next_func(real_x,time) + chol(Q)'*randn(5,1) ;
	real_data(:,i) = real_x;
	measurements(:,i) = C(real_x) + chol(R)'*randn(2,1);
end

%x and y coordinates of measurements
xx = [];
yy = [];
for i=1:outputs
	xx = [xx,measurements(1,i)*cos(measurements(2,i))];
	yy = [yy,measurements(1,i)*sin(measurements(2,i))];
end


t_start = 0;
state_count = 5;
sensor_count = 2;

%run cubature 
[estimates_cubature, covariances_cubature] = cubature(next_func,dt_between_measurements,t_start,state_count,sensor_count,...
  outputs,C,chol(Q)',chol(R)',chol(P_0)',x_0, measurements);
 
%run ukf
C = @(x,k) C(x) + 0.*k;
Q = @(k) Q;
R = @(k) R;
[estimates_ukf, covariances_ukf] = ukf(next_func, state_count, sensor_count, outputs,C,Q,R,P_0,x_0, measurements);


h = figure;
hold on;
plot(estimates_cubature(1,:),estimates_cubature(3,:),'r','LineWidth',4)
plot(estimates_ukf(1,:),estimates_ukf(3,:),'g','LineWidth',4)
plot(real_data(1,:),real_data(3,:),'k','LineWidth',2);
legend('Cubature estimate','UKF estimates','Real position')
xlabel('x')
ylabel('y')
