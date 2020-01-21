clear all;
addpath('..');
C = [1, 0.1]; %observation matrix to obtain measurement y (y = C*x + noise)
A = [-1,0.2;-0.1,-1]; %(d/dt)x(t) = A*x(t)
R_c = [0.1^2]; %continuous sensor noise
Q_c = [0.001^2,0;0,0.001^2]; %continuous process noise
P_0 = eye(2);%initial covariance
  
%next three variables dictate time step size
t_start = 0;
t_final = 20;
outputs = 10000;
 
covariances = cckf(t_start,t_final,outputs,P_0,A,C,R_c,Q_c);
	
covariances{1}%at t = 0
covariances{2000}%at t = 1.999 
format long;
covariances{end}%at t = 19.999 
