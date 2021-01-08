import pprint
import scipy
import scipy.linalg
import numpy

"""
estimate, covariance_sqrt = predictPhase(func,jacobian_func,t,P_0_sqrt,X_0,Q_root):
runs the predict portion of the dd-ekf
INPUT:
    func: x_{k+1} = f_func(x_k,t), where x_k is the state. The function's
        second argument is time t at step k ( t(k) ).
    jacobian_func: The jacobian of f_func at state x and at time (t)
        jacobian_func(x,t)
    t: current_time
    x_0: estimate of the state at time t
    P_0_sqrt: square root factor of state's covariance
        P_0 = P_0_sqrt.dot(P_0_sqrt.transpose())
    Q_root: square root of process noise covariance matrix
        Q = Q_root.dot(Q_root.transpose())
OUTPUT:
    estimate: estimate of state x after predict phase
    covariance_sqrt: square root of state's covariance matrix P
        P = covariance_sqrt.dot(covariance_sqrt.transpose())
"""
def predictPhase( func, jacobian_func, t, P_0_sqrt, x_0, Q_root ):

    x = x_0
    state_count = x.shape[0]

    estimate = func( x, t )
    jacobian = jacobian_func(x, t)

    tmp = numpy.zeros((state_count,state_count*2))
    tmp[0:state_count,0:state_count] = jacobian.dot(P_0_sqrt)
    tmp[0:state_count,state_count:] = Q_root
    Q, R = scipy.linalg.qr( tmp.transpose() )
    covariance_sqrt = R.transpose()
    covariance_sqrt = covariance_sqrt[0:state_count, 0:state_count]
    return estimate, covariance_sqrt

"""
estimate, covariance_sqrt = updatePhase(R_root,P_root,C,estimate,measurement):
runs update portion of the dd-ekf
INPUT:
    R_root: root of sensor error covariance matrix R where
        R = R_root.dot(R_root.transpose())
    P_root: root of state's covariance matrix P where
        P = P_root.dot(P_root.transpose())
    C: observation matrix
    estimate: the current estimate of the state
    measurement: sensor's measurement of the true state
OUTPUT:
    estimate: estimate of state after update phase
    covariance_sqrt: square root of state's covariance matrix P
        P = covariance_sqrt.dot(covariance_sqrt.transpose())
"""
def updatePhase( R_root, P_root, C, estimate, measurement ):
    measurement_count = C.shape[0]
    state_count = estimate.shape[0]

    tmp = numpy.zeros((state_count + measurement_count, state_count + measurement_count ))
    tmp[0:measurement_count,0:measurement_count ] = R_root
    tmp[0:measurement_count,(measurement_count):] = C.dot(P_root)
    tmp[(measurement_count):,(measurement_count):] = P_root

    Q, R = scipy.linalg.qr(tmp.transpose())

    R = R.transpose()

    X = R[ 0:measurement_count,0:measurement_count]
    Y = R[ (measurement_count):,0:(measurement_count)]
    Z = R[ (measurement_count):,(measurement_count):]

    estimate_next = estimate + Y.dot(scipy.linalg.solve(X,measurement-C.dot(estimate)))
    covariance_sqrt = Z

    return estimate_next, covariance_sqrt

"""
estimates,covariances = ddekf(func,jacobian_func,dt_between_measurements,start_time,state_count,sensor_count,measurement_count,
    C,Q_root,R_root,P_0_root,x_0,measurements):
Runs discrete-discrete Extended Kalman filter on data. The initial estimate
and covariances are at the time step before all the measurements - be
wary of the off-by-one error. If func is a linear function, then the
code is equivalent to discrete-discrete Kalman filter.

This function is meant to be used for post processing - once you have
collected all the measurements and are looking to run filtering. For real time filtering
see the implementation of this function in detail. In particular the line:

x_k_p, P_root_kp = updatePhase(R_root,P_root_km,C,x_k_m,measurements[k])

measurements[k] would be replaced with real time measurements.

INPUT:
    func: x_{k+1} = func(x_k,t) where x_k is the state. The
        function's second argument is time t (t_k) for cases when the function
        changes with time. The argument can be also used an internal 
        counter variable for func when start_time is set to zero and 
        dt_between_measurements is set to 1. 

    jacobian_func(x,t): jacobian of f_func with state x at time t

    dt_between_measurements: time distance between incoming 
        measurements. Used for incrementing time counter for each
        successive measurement with the time counter initialized with
        start_time. The time counter is fed into f_func(x,t) as t.

    start_time: the time of first measurement 

    state_count: dimension of the state

    sensor_count: dimension of observation vector

    C: observation matrix of size 'sensor_count by state_count'

    R_root: The root of sensor error covariance matrix R where
        R = R_root*(R_root'). R_root is of size 'sensor_count by 
        sensor_count'. R_root = chol(R)' is one way to derive it.

    Q_root: The root of process error covariance matrix	Q where
        Q = Q_root*(Q_root'). Q_root is of size 'state_count by 
        state_count'. Q_root = chol(Q)' is one way to derive it.

    P_0_root: The root of initial covariance matrix P_0 where
        P_0 = P_0_root*(P_root'); P_0_root is of size 'state_count by 
        state_count'. %	P_0_root = chol(P_0)' is one way to derive it.

    x_0:Initial state estimate of size 'state_count by 1'

    measurements: ith entry is ith measurement. Matrix of size 
        'sensor_count by measurement_count'
OUTPUT:
    estimates: array with 'measurement_count+1' entries where the ith
        entry the estimate of x at t_{i}
    covariances: array with 'measurements_count+1' entries where the ith
        entry is the covariance of the estimate at t_{i}
"""
def ddekf( func, jacobian_func, dt_between_measurements, start_time, state_count, sensor_count, measurement_count,
        C, Q_root, R_root, P_0_root, x_0, measurements):

    x_km1_p = x_0
    P_root_km1_p = P_0_root

    current_time = start_time

    estimates = [ x_km1_p ]
    covariances = [ P_0_root.dot(P_0_root.transpose()) ]


    for k in range(0,measurement_count):
        x_k_m, P_root_km = predictPhase(func,jacobian_func, current_time,P_root_km1_p,x_km1_p,Q_root)
        x_k_p, P_root_kp = updatePhase(R_root,P_root_km,C,x_k_m,measurements[k])

        x_km1_p = x_k_p
        P_root_km1_p = P_root_kp

        current_time = current_time + dt_between_measurements

        estimates.append(  x_km1_p  )
        covariances.append(  P_root_km1_p.dot(P_root_km1_p.transpose()) )
    return estimates, covariances
