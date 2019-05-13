function [estimate_next, covariance_sqrt] = cdekf_update_phase(R_root,P_root,C,estimate,measurement)
%runs update portion of the cd-ekf.
%[estimate_next,covariance_sqrt] = cdekf_update_phase(R_root,P_root,C,estimate,measurement)
%INPUT:
%	R_root:sensor error covariance matrix where 
%	R=R_root*(R_root')
%
%	P_root:covariance matrix from the predict phase
%	P_predict = P_root*(root')
%
%	C: observation matrix
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

	measurement_count = size(C,1); 
	state_count = size(C,2);
	
	tmp = zeros(state_count+measurement_count, state_count+measurement_count);
	tmp(1:measurement_count,1:measurement_count) = R_root;
	tmp(1:measurement_count,(measurement_count+1):end) = C*P_root;
	tmp((1+measurement_count):end,(1+measurement_count):end) = P_root;

	[Q,R] = qr(tmp');
	R = R';

	X = R(1:measurement_count,1:measurement_count); 
	Y = R((1+measurement_count):end,(1:measurement_count));
	Z = R((1+measurement_count):end,(1+measurement_count):end);	
	
	%Kalman_gain = Y*inv(X);
	estimate_next = estimate + Y*(X\(measurement - C*estimate));
	covariance_sqrt = Z;
end