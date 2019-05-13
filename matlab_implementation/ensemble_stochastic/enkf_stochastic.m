function [estimates,covariances] = enkf_stochastic(func,state_count,sensor_count,measurement_count,ensemble_size,C,R,Q,P_0,x_0, measurements)
%Runs discrete-discrete Extended Kalman filter on data. The initial
%estimate and covariances are at the time step before all the 
%measurements - be wary of the off-by-one error. If f_func is a
%linear function the code is equivalent to discrete-discrete Kalman
%filter. 
%[estimates, covariances] = ddekf(f_func,jacobian_func,state_count,sensor_count,measurement_count,C,Q,R_root,P_0_root,x_0, measurements)
%INPUT:
%	f_func: x_{k+1} = f_func(x_k) where x_k is the state. 
%
%	state_count: dimension of the state
%
%	sensor_count: dimension of observation vector
%
%	measurement_count: number of measurements
%
%	ensemble_size: Ensemble size of the ensemble filter
%
%	C: observation matrix of size 'sensor_count by state_count'
%
%	R: Sensor error covariance
%
%	Q: Process noise covariance
%
%	P: Initial covariance
%
%	x_0: Initial state estimate of size 'state_count by 1'
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


assert(size(C,1) == sensor_count & size(C,2) == state_count); 
assert(size(R,1) == sensor_count & size(R,2) == sensor_count);
assert(size(P_0,1) == state_count & size(P_0,2) == state_count);
assert(size(x_0,1) == state_count & size(x_0,2) == 1); 
assert(size(measurements,1) == sensor_count & size(measurements,2) == measurement_count);
assert(size(func(zeros(state_count,1)),1) == state_count & size(func(zeros(state_count,1)),2) == 1);

covariances = cell(1,measurement_count);
estimates = zeros(state_count,measurement_count);

ensemble = mvnrnd(zeros(state_count,1),P_0,ensemble_size)' + x_0;

for ii=1:measurement_count
	%propagate
	propagated = func(ensemble);
	sampled_process_noise = mvnrnd(zeros(state_count,1),Q,ensemble_size)';
	propagated = propagated + sampled_process_noise;
	x_mean = (1/ensemble_size).*sum(propagated,2);
	propagated_covariance = (propagated - x_mean);
	propagated_covariance = (1/(ensemble_size-1)).*propagated_covariance*propagated_covariance';
	Kalman_gain = propagated_covariance*C'/(C*propagated_covariance*C' + R);
	
	%draw a stastistically consistent observation set:
	u_samples = mvnrnd(zeros(sensor_count,1),R,ensemble_size)';
	y_propagated = C*propagated; 
	
	updated_ensemble = propagated  + Kalman_gain*(measurements(:,ii) + u_samples - y_propagated);
	ensemble = updated_ensemble;
	
	estimates(:,ii) = (1/ensemble_size).*sum(ensemble,2);
	covariances{ii} = (1/(ensemble_size-1)).*((ensemble-estimates(:,ii)))*((ensemble-estimates(:,ii)))';
end

end
