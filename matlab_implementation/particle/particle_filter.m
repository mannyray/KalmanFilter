function [estimates, covariances] = particle_filter(f_func,jacobian_func,dt_between_measurements,start_time,state_count,sensor_count,measurement_count,particle_count,particle,C,Q_root,R_root,P_0_root,x_0, measurements)
%Particle filter using Extended Kalman filter for local linearization
%and Sequential Importance Resampling (SIR) for the resampling.
%estimate and covariances are at the time step before all the
%measurements - be wary of the off-by-one error. If f_func is a
%linear function the code is equivalent to discrete-discrete Kalman
%filter.
%[estimates, covariances] = particle_filter(f_func,jacobian_func,state_count,sensor_count,measurement_count,particle_count,particle,C,Q_root,R_root,P_0_root,x_0, measurements)
%INPUT:
%       f_func: x_{k+1} = f_func(x_k,t) where x_k is the state. The
%       function's second argument is time t for cases when the function
%       changes with time. The argument can be also used an internal
%       counter variable for f_func when start_time is set to zero and
%       dt_between_measurements is set to 1.
%
%       jacobian_func(x,t): jacobian of f_func with state x at time t
%
%       dt_between_measurements: time distance between incoming
%       measurements. Used for incrementing time counter for each
%       successive measurement with the time counter initialized with
%       start_time. The time counter is fed into f_func(x,t) as t.
%	(usually set to 1)
%
%       start_time: the time of first measurement
%
%       state_count: dimension of the state
%
%       sensor_count: dimension of observation vector
%
%	measurement_count: The total amount of measurements.
%
%	particle_count: The amount of particles in the particle_filter.
%
%	particle: `state_count X particle_count` matrix where the i^th
%	column is the i^th particle. Can be initialized as
%	`particle = mvnrnd(x_0,P_0,particle_count)';` on user end. The 
%	particle weights will be initialized as uniform.
%
%       C: observation matrix of size 'sensor_count by state_count'
%
%       R_root: The root of sensor error covariance matrix R where
%       R = R_root*(R_root'). R_root is of size 'sensor_count by
%       sensor_count'. R_root = chol(R)' is one way to derive it.
%
%       Q_root: The root of process error covariance matrix     Q where
%       Q = Q_root*(Q_root'). Q_root is of size 'state_count by
%       state_count'. Q_root = chol(Q)' is one way to derive it.
%
%       P_0_root: The root of initial covariance matrix P_0 where
%       P_0 = P_0_root*(P_root'); P_0_root is of size 'state_count by
%       state_count'. % P_0_root = chol(P_0)' is one way to derive it.
%	Particles will be initialized with P_0 as the initial covariance.
%
%       x_0:Initial state estimate of size 'state_count by 1'
%
%       measurements: ith column is ith measurement. Matrix of size
%       'sensor_count by measurement_count'
%
%OUTPUT:
%       estimates: 'state_count by measurement_count+1'
%       ith column is ith estimate. first column is x_0
%
%       covariances: cell of size 'measurement_count+1' by 1
%       where each entry is the P covariance matrix at that time
%       Time is computed based on dt_between_measurements

%basic input checking
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
assert(size(particle,1)==state_count&&...
	size(particle,2)==particle_count);


estimates = zeros(state_count,measurement_count+1);
covariances = cell(measurement_count+1,1);
particle_covariance = cell(particle_count,1);
%covariance is set to be the same for all initial particles
for particle_index=1:particle_count
	particle_covariance{particle_index} = P_0_root*P_0_root';
end
%initialized as uniform
particle_weight = (1/particle_count).*ones(particle_count,1);

estimates(:,1) = x_0;
covariances{1} = P_0_root*P_0_root';

for k=1:measurement_count
	%project EKF in each particle i
	%this (of course) can be parallelized as other for loops here can too
	for particle_index=1:particle_count
		%needed just for the covariance - estimate still computed within ddekf_predict_phase - can be removed to further optimize
		[covariance_sqrt] = predict_phase(f_func,jacobian_func,k,particle_covariance{particle_index},particle(:,particle_index),Q_root);

		particle(:,particle_index) = f_func(particle(:,particle_index),k) + mvnrnd(zeros(state_count,1),Q_root*Q_root')';
		particle_covariance{particle_index} = covariance_sqrt;
	end

	%update EKF in each particle i
	for particle_index=1:particle_count
		[covariance_sqrt] = update_phase(R_root,particle_covariance{particle_index},C,particle(:,particle_index),measurements(k));	

		particle(:,particle_index) = estimate;
		particle_covariances{particle_index} = covariance_sqrt;
	end

	%update particle weights
	for weight_index=1:particle_count
		particle_weight(weight_index) = particle_weight(weight_index)*normpdf(measurements(k),C*particle(:,weight_index),R_root*R_root');
	end
	%https://stats.stackexchange.com/questions/201545/likelihood-calculation-in-particle-filtering
	%normalize particle weights
	particle_weight = (1/sum(particle_weight)).*particle_weight;

	%compute and store mean estimate and mean covariance
	mean_estimate = zeros(state_count,1);
	covariance = zeros(state_count,state_count);
	for particle_index=1:particle_count
		mean_estimate = mean_estimate + particle_weight(particle_index)*particle(:,particle_index);
	end
	for particle_index=1:particle_count
		covariance = covariance + particle_weight(particle_index).*( (particle(:,particle_index) - mean_estimate)*(particle(:,particle_index) - mean_estimate)'); 
	end
	estimates(:,k+1) = mean_estimate;
	covariances{k+1} = covariance;	

	%resample
	%compute cumulitive profile using a rolling sum
	%(cumulitive distribution function [cdf] - increasing
	%function with 0 at -infinity and 1 at +infinity).
	cumulitive_profile = zeros(particle_count,1);
	cumulitive_profile(1) = particle_weight(1);
	for weight_index=2:particle_count
		cumulitive_profile(weight_index) = cumulitive_profile(weight_index-1) + particle_weight(weight_index); 
	end
	cumulitive_profile = (1/cumulitive_profile(end)).*cumulitive_profile;%normalize

	%use uniform distribution(rand) to select from [cdf]
	%cumulitive_profile to resample/update the particles and their weights.
	%The particle with the biggest weight will have the biggest y-axis
	%'gain' in the cdf plot and therefore chosen_index will more frequently
	%select this particle and produce more particles in particle_new
	%with uniform particle weight (1/particle_count).
	particle_new = zeros(state_count,particle_count); 
	particle_covariance_new = cell(state_count,1);
	for weight_index=1:particle_count
		chosen_index = (find( cumulitive_profile > rand,1) - 1);
		particle_new(:,weight_index) = particle(:,weight_index);
		particle_covariance_new{weight_index} = particle_covariance{weight_index};
	end
	particle_weight = (1/particle_count).*ones(particle_count,1);
	particle = particle_new;
	particle_covariance = particle_covariance_new;
end

end
