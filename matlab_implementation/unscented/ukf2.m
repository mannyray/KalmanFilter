function [estimates, covariances] = ukf2(f_func, state_count, sensor_count, measurement_count,C,Q,R,P_0,x_0, measurements)
%_Faster_ Unscented Kalman filter implementation.
%[estimates, covariances] = ukf2(f_func, state_count, sensor_count, measurement_count,C,Q,R,P_0,x_0, measurements)
%INPUT
%	f_func: function when given a state, produces the next state 
%	x_{k+1} = f_func(x_{k})
%	
%	state_count: length of state
%
%	sensor_count:amount of measurements
%
%	measurement_count: total number of measurements
%
%	C: observation function that takes in state X_k and index k and returns 
%	appropriate measurement: observation = C(X_k,k) 
%
%	Q: function that takes in argument index k and returns state_count by
%	state_count process noise matrix: Q_k = Q(k)
%
%	R: function that takes in argument index k and returns sensor_count by
%	sensor_count sensor noise covariance: R_k = R(k)
%
%	P_0: state_count by state_count initial covariance
%
%	x_0: initial estimate of state (state_count by 1)
%
%	measurements: sensor_count by measurement_count matrix where the k^th
%	column is the k^th measurement
%
%OUTPUT
%	estimates: state_count by (measurement_count + 1) matrix where the 
%	k^th column is the k^th estimate (first column is x_0)
%
%	covariances: 1 by (measurement_count + 1) cell count where the kth 
%	index is the computed (state_count by state_count) covariance matrix
%	(first entry is P_0)


assert(state_count == size(x_0,1));
for i=1:measurement_count
	assert(size(C,1) == sensor_count & size(C,2) == state_count);
	assert(size(Q,1) == state_count & size(Q,2) == state_count);
	assert(size(R,1) == sensor_count & size(R,2) == sensor_count);
end
assert(size(P_0,1) == state_count & size(P_0,2) == state_count);
assert(size(x_0,1) == state_count & size(x_0,2) == 1); 
assert(size(measurements,1) == sensor_count & size(measurements,2) == measurement_count);
assert(size(f_func(zeros(state_count,1)),1) == state_count & size(f_func(zeros(state_count,1)),2) == 1);


P_update = P_0;
x_update = x_0;
covariances = cell(1,measurement_count+1);
estimates = zeros(state_count,measurement_count+1);
covariances{1,1} = P_0;
estimates(:,1) = x_0;

sigma_points = zeros(state_count,state_count*2+1);
propagated = zeros(state_count,state_count*2,+1);
measurement_points = zeros(sensor_count,state_count*2+1);
cross_covariance = zeros(state_count,sensor_count);
Z_cov = zeros(sensor_count,sensor_count);
P_cov = zeros(state_count,state_count);

alpha = 1;
beta = 2;
KK = 0;
lambda = alpha^2*(state_count+KK) - state_count;

n_0_m = (lambda/(state_count + lambda));
n_i_m = 0.5*(1/(state_count + lambda));
n_0_c = (lambda/(state_count+lambda) + 1 - alpha^2 + beta); 
n_i_c = n_i_m;

for ii=1:measurement_count
	P_sqrt = chol(P_update)';
	sigma_points = x_update +  (sqrt(state_count + lambda).*[zeros(state_count,1),P_sqrt, -P_sqrt]);

	propagated = f_func(sigma_points);
	measurement_points = C*sigma_points;

	x_predict = n_0_m.*propagated(:,1) + n_i_m.*(sum(propagated(:,2:end),2));
	z_predict = n_0_m.*measurement_points(:,1) + n_i_m.*(sum(measurement_points(:,2:end),2));
    

	Z_cov = n_0_c.*(z_predict - measurement_points(:,1))*(z_predict - measurement_points(:,1))';
	P_cov = n_0_c.*(x_predict - propagated(:,1))*(x_predict-propagated(:,1))';
	cross_covariance = n_0_c.*(x_predict - propagated(:,1))*(z_predict - measurement_points(:,1))';

    
	diff_z = measurement_points(:,2:end) - z_predict;
	diff_x = propagated(:,2:end) - x_predict;
	Z_cov = Z_cov + n_i_c.*diff_z*(diff_z');
	P_cov = P_cov + n_i_c.*diff_x*(diff_x');
    	cross_covariance = cross_covariance + n_i_c.*(diff_x*(diff_z'));

	P_cov = P_cov + Q;
	P_cov =  0.5.*(P_cov + P_cov');
	Z_cov = Z_cov + R;
	Z_cov = 0.5*(Z_cov + Z_cov');

	%update phase
	Kalman_gain = cross_covariance*inv(Z_cov);
	x_update = x_predict + Kalman_gain*(measurements(:,ii) - z_predict);

	P_update = P_cov - Kalman_gain*Z_cov*(Kalman_gain');
	covariances{1,ii+1} = P_update;
	estimates(:,ii+1) = x_update;
end
