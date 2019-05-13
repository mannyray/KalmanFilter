function [estimates, covariances ] = cdekf(f_func,jacobian_func,dt_between_measurements,rk4_steps, start_time, state_count, sensor_count, measurement_count,C,Q_root,R_root,P_0_root,x_0, measurements)
%Runs continuous-discrete Extended Kalman filter on data. The initial
%estimate and covariances are at the time step before all the 
%measurements - be wary of the off-by-one error. If f_func is a
%linear function the code is equivalent to continuous-discrete Kalman
%filter. 
%[estimates, covariances] = cdekf(f_func,jacobian_func,dt_between_measurements,rk4_steps, start_time, state_count, sensor_count, measurement_count,C,Q_root,R_root,P_0_root,x_0, measurements)
%INPUT:
%	f_func: \dot{x} = f_func(x,t) where x is the state
% 	and lhs is the derivative of x at time t. The derivative is a 
%	function of both the time and state.
%
%	jacobian_func(x,t): jacobian of f_func with state x at time t
%
%	rk4_steps: amount of steps taken in predict phase in rk4 scheme
%
%	dt_between_measurements: time distance between incoming
%	measurements
%
%	start_time: the time of first measurement 
%
%	state_count: dimension of the state
%
%	sensor_count: dimension of observation matrix
%
%	C: observation matrix of size 'sensor_count by state_count'
%
%	R_root: The root of sensor error covariance matrix R where
%	R = R_root*(R_root'). R_root is of size 'sensor_count by 
%	sensor_count'. R_root = chol(R)' is one way to derive it.
%
%	Q_root: The root of process error covariance matrix Q where
%	Q = Q_root*(Q_root'). Q_root is of size 'state_count by 
%	state_count'. Q_root = chol(Q)' is one way to derive it.
%
%	P_0_root: The root of initial covariance matrix P_0 where
%	P_0 = P_0_root*(P_root'); P_0_root is of size 'state_count by 
%	state_count'. P_0_root = chol(P_0)' is one way to derive it.
%
%	x_0:Initial state estimate of size 'state_count by 1'
%
%	measurements: ith column is ith measurement. Matrix of size 
%	'sensor_count by measurement_count'
%
%OUTPUT:
%	estimates: 'state_count'X'measurement_count+1'
%	ith column is ith estimate. first column is x_0
%
%	covariances: cell of size 'measurement_count+1' by 1
%	where each entry is the P covariance matrix at that time.
%	Time is computed based on dt_between_measurements

	%make sure input is valid
	assert(size(P_0_root,1)==state_count &&...
		size(P_0_root,2)==state_count);
	assert(size(C,1)==sensor_count && size(C,2)==state_count);
	assert(size(Q_root,1)==state_count &&...
		size(Q_root,2)==state_count);
	assert(size(x_0,1)==state_count && size(x_0,2)==1);
	assert(size(R_root,1)==sensor_count &&...
		size(R_root,2)==sensor_count);
	test = f_func(x_0,start_time);
	assert(size(test,1)==state_count && size(test,2)==1);
	test = jacobian_func(x_0,start_time);
	assert(size(test,1)==state_count && size(test,2)==state_count);
	assert(size(measurements,1)==sensor_count &&...
		size(measurements,2)==measurement_count);

	current_time = start_time;
	h = dt_between_measurements/rk4_steps;
	current_time = 0;
	rk4_step = dt_between_measurements;

	x_km1_p = x_0;
	P_root_km1_p = P_0_root;

	estimates = zeros(state_count,measurement_count+1); 
	covariances = cell(measurement_count + 1, 1);
	
	estimates(:,1) = x_0; 
	covariances{1,1} = P_0_root*P_0_root';
	
	for k=1:measurement_count
		%loop through all measurements and follow through predict 
		%and update phase
		[x_k_m,P_root_km] = cdekf_predict_phase(f_func,...
			jacobian_func,h,dt_between_measurements,current_time,...
			P_root_km1_p,x_km1_p,Q_root);
		[x_k_p,P_root_kp] = cdekf_update_phase(R_root,...
			P_root_km,C,x_k_m,measurements(:,k));

		x_km1_p = x_k_p;
		P_root_km1_p = P_root_kp;

		%store results
		estimates(:,k+1) = x_k_p;
		covariances{k+1,1} = P_root_kp*(P_root_kp');

		current_time = current_time + dt_between_measurements; 
	end
end
