# Steady-state Kalman Filter

The steady state of continuous\_continuous Kalman filter which is equivalent to when time derivative of `P(t)` goes to zero and reduces the Differential Riccati equation to the Continuous Algebraic Riccati equation (works in Matlab only). The code example is available in `matlab_implementation/steady_state/examples/ss_example.m` and the code for the filter is in `matlab_implementation/steady_state` of [https://github.com/mannyray/KalmanFilter](https://github.com/mannyray/KalmanFilter):

```
clear all;
addpath('..');
C = [1, 0.1]; %observation matrix to obtain measurement y (y = C*x + noise)
A = [-1,0.2;-0.1,-1]; %(d/dt)x(t) = A*x(t)
R_c = [0.1^2]; %continuous sensor noise
Q_c = [0.001^2,0;0,0.001^2]; %continuous process noise

P_steady = sskf(A,C,Q_c,R_c)
```

produces

```
ans = 

1.0e-06 *

 0.504888842491879   0.024508560851831
 0.024508560851831   0.497548868161834
```
