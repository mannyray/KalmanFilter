clear all;

addpath('../continuous_discrete');
addpath('../discrete_discrete');

P_cell = {[5,0.5;0.5,3], [1 2 6;2 9 0.5; 6,0.5,0.1]};
for i=1:length(P_cell)
	P = P_cell{1,i};
	state_count = size(P,1);
	r = 0.001;
	P_root = chol(P)';
	assert(max(max(abs(P_root*P_root'-P)))<0.0000000001);
	R_root = sqrt(r);
	measurement = 122.01;
	estimate = randn(state_count,1);
	C = randn(1,state_count); 

	[estimate_next, covariance_sqrt] = ddekf_update_phase(R_root,P_root,C,estimate,measurement);
	[estimate_next2, covariance_sqrt2] = cdekf_update_phase(R_root,P_root,C,estimate,measurement);

	%manual test (the classic udpate filter algorithm without square root algorithm)
	K_k = (P_root*P_root')*(C')*inv(C*(P_root*P_root')*C'+R_root*R_root);
	coV = (eye(state_count) - K_k*C)*(P_root*P_root');
	estimate_new = estimate + K_k*(measurement - C*estimate);

	assert(norm(estimate_new - estimate_next) < 0.000000000001);
	assert(max(max(abs(covariance_sqrt*covariance_sqrt' - coV))) < 0.0000000001)

end
disp('update phase tested succesfully')
