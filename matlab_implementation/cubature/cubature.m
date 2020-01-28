function [estimates, covariances ] = cubature(f_func,dt_between_measurements,start_time,state_count,sensor_count,measurement_count,C_func,Q_root,R_root,P_0_root,x_0, measurements)
%Runs Cubature Kalman filter on data. The initial
%estimate and covariances are at the time step before all the 
%measurements - be wary of the off-by-one error. If f_func is a
%linear function the code is equivalent to discrete-discrete Kalman
%filter. 
%[estimates, covariances] = cubature(f_func,jacobian_func,state_count,sensor_count,measurement_count,C,Q_root,R_root,P_0_root,x_0, measurements)
%INPUT:
%	f_func: x_{k+1} = f_func(x_k,t) where x_k is the state. The
%	function's second argument is time t for cases when the function
%	changes with time. The argument can be also used an internal 
%	counter variable for f_func when start_time is set to zero and 
%	dt_between_measurements is set to 1. 
%
%	dt_between_measurements: time distance between incoming 
%	measurements. Used for incrementing time counter for each
%	successive measurement with the time counter initialized with
%	start_time. The time counter is fed into f_func(x,t) as t.
%
%	start_time: the time of first measurement 
%
%	state_count: dimension of the state
%
%	sensor_count: dimension of observation vector
%
%	C_func: observation matrix function that converts single argument vector 
%	of dimension of 'state_count' by 1 to system measurement of size 'sensor_count' by 1
%
%	R_root: The root of sensor error covariance matrix R where
%	R = R_root*(R_root'). R_root is of size 'sensor_count by 
%	sensor_count'. R_root = chol(R)' is one way to derive it.
%
%	Q_root: The root of process error covariance matrix	Q where
%	Q = Q_root*(Q_root'). Q_root is of size 'state_count by 
%	state_count'. Q_root = chol(Q)' is one way to derive it.
%
%	P_0_root: The root of initial covariance matrix P_0 where
%	P_0 = P_0_root*(P_root'); P_0_root is of size 'state_count by 
%	state_count'. %	P_0_root = chol(P_0)' is one way to derive it.
%
%	x_0:Initial state estimate of size 'state_count by 1'
%
%	measurements: ith column is ith measurement. Matrix of size 
%	'sensor_count by measurement_count'
%
%OUTPUT:
%	estimates: 'state_count by measurement_count+1'
%	ith column is ith estimate. first column is x_0
%
%	covariances: cell of size 'measurement_count+1' by 1
%	where each entry is the P covariance matrix at that time
%	Time is computed based on dt_between_measurements

	%make sure input is valid
	assert(size(P_0_root,1)==state_count &&...
	 size(P_0_root,2)==state_count);
	%assert(size(C,1)==sensor_count && size(C,2)==state_count);
	assert(size(Q_root,1)==state_count &&...
	 size(Q_root,2)==state_count);
	assert(size(x_0,1)==state_count && size(x_0,2)==1);
	assert(size(R_root,1)==sensor_count &&...
	 size(R_root,2)==sensor_count);
	test = f_func(x_0,start_time);
	assert(size(test,1)==state_count && size(test,2)==1);
	%test = jacobian_func(x_0,start_time);
	%assert(size(test,1)==state_count && size(test,2)==state_count);
	assert(size(measurements,1)==sensor_count &&...
	 size(measurements,2)==measurement_count);

	x_km1_p = x_0;
	P_root_km1_p = P_0_root;

	estimates = zeros(state_count,measurement_count + 1); 
	covariances = cell(measurement_count + 1, 1);
	
	estimates(1:state_count,1) = x_0; 
	covariances{1,1} = P_0_root*P_0_root';

	current_time = start_time;	
	
	for k=1:measurement_count
		%loop through all measurements and follow through predict 
		%and update phase
		[x_k_m,P_root_km] = cubature_predict_phase(f_func,...
			current_time,P_root_km1_p,x_km1_p,...
			Q_root);
   		[x_k_p,P_root_kp] = cubature_update_phase(R_root,...
   			P_root_km,C_func,x_k_m,measurements(:,k));
		
		x_km1_p = x_k_p;
		P_root_km1_p = P_root_kp;

		%store results
		estimates(:,k+1) = x_k_p;
		covariances{k+1,1} = P_root_kp*(P_root_kp');

		current_time = current_time + dt_between_measurements; 
	end
end
