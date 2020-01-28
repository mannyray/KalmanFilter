function [estimate_next, covariance_sqrt] = cubature_update_phase(R_root,P_root,C_func,estimate,measurement)
%runs update portion of the dd-ekf.
%[estimate_next,covariance_sqrt] = cubature_update_phase(R_root,P_root,C,estimate,measurement)
%INPUT:
%	R_root:sensor error covariance matrix where 
%	R=R_root*(R_root')
%
%	P_root:covariance matrix from the predict phase
%	P_predict = P_root*(root')
%
%	C_func: observation matrix function that converts single argument vector 
%	of dimension of 'state_count' by 1 to system measurement of size 'sensor_count' by 1
%
%	estimate: estimate from predict phase
%
%	measurement: measurement picked up by sensor
%
%OUTPUT:
%	estimate_next: estimate incorporating measurement 
%
%	covariance_sqrt: square root of covariance at the 
%	end of update phase. Actual covariance is 
%	covariance_sqrt*(covariance_sqrt')

	state_count = length(estimate);
	sensor_count = length(R_root);
	
	sigma_points = zeros(state_count,state_count*2);
	sigma_points(:,1:state_count) = estimate + P_root*sqrt(state_count).*eye(state_count,state_count);
	sigma_points(:,(state_count+1):end) = estimate - P_root*sqrt(state_count).*eye(state_count,state_count);

	measurement_points = zeros(sensor_count,state_count*2);
	for i=1:(2*state_count)
		measurement_points(:,i) = C_func(sigma_points(:,i));
	end
	average_measurement = (1/(2*state_count))*sum(measurement_points,2);
	
	Z_cov = 0;
	cross_covariance = 0;
	for i=1:(state_count*2)
		Z_cov = (average_measurement - measurement_points(:,i))*(average_measurement - measurement_points(:,i))';
		cross_covariance = (estimate  - sigma_points(:,i))*(average_measurement - measurement_points(:,i))';
	end
	Z_cov = Z_cov + R_root*R_root';
	Kalman_gain = cross_covariance*inv(Z_cov);

	estimate_next = estimate + Kalman_gain*(measurement - average_measurement);
	covariance = P_root*P_root' - Kalman_gain*Z_cov*(Kalman_gain');
	covariance_sqrt = chol(covariance)';
end
