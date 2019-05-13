addpath('../continuous_discrete');

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

	[estimate, covariance_sqrt] = cdekf_predict_phase(f_func,jacobian_func,h,del_t,start_time,P_0_sqrt,X_0,Q_root);

	cur_time = 0;
	p = P;
	x = X_0;
	phi = eye(state_count);
	phi_sum = zeros(state_count);
	while(cur_time < del_t)
		phi_next = phi + (h).*(jacobian_func(x,cur_time)*phi); 
		phi_sum = phi_sum + (h/2).*( phi*Q*phi' + phi_next*Q*phi_next');

		x = x + h*f_func(x,cur_time);
		p = p + h.*( jacobian_func(x,cur_time)*p + p*jacobian_func(x,cur_time)' + Q_root*Q_root);
		cur_time = cur_time + h;
	
		phi = phi_next;
	end
	phi_sum = phi_sum + phi*P*phi';
	
	assert(max(max(abs(covariance_sqrt*(covariance_sqrt')-p)))<0.000000001);
	
end
