# Continuous-continuous Kalman Filter

The states evolves continuously in time and the measurements are also continuous in time. In this implementation only the covariance `P(t)` is solved for through the Differential Riccati equation and that only works for linear systems. A basic numerical scheme is used so the code should only be used for small state count and care taken in setting timestep. The code example is available in `matlab_implementation/continuous_continuous/examples/cc_example.m` and the code for the filter is in `matlab_implementation/continuous_continuous` of [https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter). 

```
clear all;
addpath('../');
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
```

produces

```
ans =

   8.2441e-05   1.1739e-04
   1.1739e-04   1.8665e-04

ans =

   5.04888842493018e-07   2.45085608502966e-08
   2.45085608502966e-08   4.97548868164108e-07
```
