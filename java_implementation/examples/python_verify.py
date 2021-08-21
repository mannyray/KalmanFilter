import sys
sys.path.insert(0,'../../python_implementation/discrete_discrete')
from ddekf import *
import numpy as np

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

measurements_file = open('measurements_java.txt')
measurement_lines = measurements_file.readlines()

measurements = []
for line in measurement_lines:
    measurement = np.array([[float(line)]])
    measurements.append(measurement)


estimates, covariances = ddekf( func, jacobian_func, dt, t_start, 
     state_count, sensor_count, measurement_count, C, Q_root, R_root, P_0_sqrt, x_0, measurements)

#load in what the java code computed
estimates_file = open('estimates_java.txt')
estimates_lines = estimates_file.readlines()

estimates_java = []
for line in estimates_lines:
    arr = np.fromstring(line, dtype=float, sep=',')
    estimate = np.zeros((2,1))
    estimate[0,0] = arr[0]
    estimate[1,0] = arr[1]
    estimates_java.append(estimate)

covariance_file = open('covariances_java.txt')
covariance_lines = covariance_file.readlines()

covariances_java = []
for line in covariance_lines:
    arr = np.fromstring(line, dtype=float, sep=',')
    covariance = np.zeros((2,2))
    covariance[0,0] = arr[0]
    covariance[1,0] = arr[1]
    covariance[0,1] = arr[2]
    covariance[1,1] = arr[3]
    covariances_java.append(covariance)

print("difference between Java and Python ddekf implementation:")
max_relative_error_estimates = 0
for index in range(1,len(estimates_java)):
    estimate_python = estimates[index]
    estimate_java = estimates_java[index]
    max_relative_error_estimates = max(max_relative_error_estimates,(np.linalg.norm(estimate_python-estimate_java)/np.linalg.norm(estimate_python)))
print("max relative error for estimates: " + str(max_relative_error_estimates))

max_relative_error_covariances = 0
for index in range(1,len(covariances_java)):
    covariance_python = covariances[index]
    covariance_java = covariances_java[index]
    max_relative_error_covariances = max(max_relative_error_covariances,(np.linalg.norm(covariance_python-covariance_java)/np.linalg.norm(covariance_python)))
print("max relative error for covariances: " + str(max_relative_error_covariances))

