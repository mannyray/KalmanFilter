import sys
sys.path.insert(0,'..')
from ddekf import *
import numpy as np


np.random.seed(1)


C = np.array([[1, 0.1]])
A = np.array([[-1,0.2],[-0.1,-1]])
R_c = np.array( [[0.1*0.1]] )
Q_c = np.array( [[0.001*0.001,0],[0,0.001*0.001]])
P_0 = np.array([[1,0],[0,1]])
x_0 = np.array([[100,80]]).transpose()
state_count = 2
sensor_count = 1

t_start = 0
t_final = 20
measurement_count = 10000
dt = (t_final-t_start)/measurement_count

A_d = np.array([[1,0],[0,1]]) + dt*A
Q_d = dt*Q_c;
R_d = (1/dt)*R_c;

Q_root = np.linalg.cholesky(Q_d).transpose()
R_root = np.linalg.cholesky(R_d).transpose()
P_0_sqrt = np.linalg.cholesky(P_0).transpose()

def func(x,t):
    return A_d.dot(x)

def jacobian_func(x,t):
    return A_d

ideal_data = []
process_noise_data = []
measurements = []
x = x_0
x_noise = x_0
for ii in range(0,measurement_count):
    x = func(x,ii)
    ideal_data.append( x )
    x_noise = func(x_noise,ii) + Q_root.transpose().dot(np.random.randn(state_count,1))
    process_noise_data.append(x_noise)
    measurements.append(C.dot(x_noise) + R_root.transpose().dot(np.random.randn(sensor_count,1)))

estimates, covariances = ddekf( func, jacobian_func, dt, t_start,
    state_count, sensor_count, measurement_count, C, Q_root, R_root, P_0_sqrt, x_0, measurements)

print( covariances[len(covariances)-1])

