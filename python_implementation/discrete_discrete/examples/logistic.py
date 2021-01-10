import sys
sys.path.insert(0,'..')
from ddekf import *
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(1)

rate = 0.01
max_pop = 100.0
def func(x,t):
    return x + rate*x*(1 -(1/max_pop)*x)

def jacobian_func(x,t):
    return 1 + rate - (2*rate/max_pop)*x


x_0 = np.zeros((1,1))
x_0[0][0] = max_pop/2
P_0 = np.array([[1]])
Q = np.array([[1]])
R = np.array([[3]])
C = np.array([[1]])

Q_root = np.linalg.cholesky(Q).transpose()
R_root = np.linalg.cholesky(R).transpose()
P_0_sqrt = np.linalg.cholesky(P_0).transpose()
Q_root = np.array([[1]])

start_time = 0
finish_time = 10
sensor_count = 1
state_count = 1
measurement_count = 1000
dt_between_measurements = (finish_time - start_time)/measurement_count

x = x_0
x_noise = x_0
ideal_data = []
process_noise_data = []
measurements = []
times = []
for ii in range(0,measurement_count):
    current_time = ii*dt_between_measurements
    times.append(current_time)
    x = func(x,current_time)
    ideal_data.append( x )
    x_noise = func(x_noise,current_time) + Q_root.transpose().dot(np.random.randn(state_count,1))
    process_noise_data.append(x_noise)
    measurements.append(C.dot(x_noise) + R_root.transpose().dot(np.random.randn(sensor_count,1)))

estimates, covariances = ddekf( func, jacobian_func, dt_between_measurements, start_time, state_count, sensor_count, measurement_count, C, Q_root, R_root, P_0_sqrt, x_0, measurements)

measurements_flat = [ x[0][0] for x in measurements ]
estimates_flat = [ x[0][0] for x in estimates ]
process_noise_data_flat = [ x[0][0] for x in process_noise_data ]

line1, = plt.plot(times,measurements_flat,color='blue',label='Measurement')
line2, = plt.plot(times,process_noise_data_flat,color='red',label='Real data')
line3, = plt.plot(times,estimates_flat[1:],color='orange',label='Estimate')
plt.legend(handles=[line1,line2,line3])
plt.ylabel('Population')
plt.xlabel('Time')
plt.plot()
plt.savefig('logistic2.png')
plt.show()
plt.close()

line1, = plt.plot(times[0:99],measurements_flat[0:99],color='blue',label='Measurement')
line2, = plt.plot(times[0:99],process_noise_data_flat[0:99],color='red',label='Real data')
line3, = plt.plot(times[0:99],estimates_flat[1:100],color='orange',label='Estimate')
plt.legend(handles=[line1,line2,line3])
plt.ylabel('Population')
plt.xlabel('Time')
plt.plot()
plt.savefig('logistic1.png')
plt.show()
plt.close()

