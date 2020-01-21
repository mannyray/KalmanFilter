clear all;
addpath('..');
C = [1, 0.1]; %observation matrix to obtain measurement y (y = C*x + noise)
A = [-1,0.2;-0.1,-1]; %(d/dt)x(t) = A*x(t)
R_c = [0.1^2]; %continuous sensor noise
Q_c = [0.001^2,0;0,0.001^2]; %continuous process noise

P_steady = sskf(A,C,Q_c,R_c)
