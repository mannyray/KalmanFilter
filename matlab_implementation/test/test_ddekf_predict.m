addpath('../discrete_discrete');

P_cell = {[5,0.5;0.5,3], [1 2 6;2 9 0.5; 6,0.5,0.1]};
for i=1:length(P_cell)
	P = P_cell{1,i}:
	state_count = size(P,1);
	r = 0.001;
	P_0_sqrt = chol(P)';
	X_0 = randn(state_count,1);
	Q_root = 1*eye(state_count);
	Q = Q_root*Q_root;
	assert(max(max(abs(P_0_sqrt*P_0_sqrt'-P)))<0.000000001);

	A = randn(state_count,state_count);
	f_func = @(x,t) A*x + 0*t;
	jacobian_func = @(x,t) A + 0.*x + 0*t;
	
	h = 0.0001;
	del_t = 2;
	start_time = 0;

	[estimate, covariance_sqrt] = ddekf_update_phase(R_root,P_root,C,estimate,measurement);

	cur_time = 0;
	p = A*P*(A') + Q;
	
	assert(max(max(abs(covariance_sqrt*(covariance_sqrt')-p)))<0.000000001);
	
end
