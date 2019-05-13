function [estimate, covariance_sqrt] = ddekf_predict_phase(f_func,jacobian_func,t,P_0_sqrt,x_0,Q_root)
%runs predict portion of update of the dd-ekf
%[estimate, covariance_sqrt] ddekf_predict_phase(f_func,jacobian_func,t,P_0_sqrt,x_0,Q_root)
%INPUT:
%	f_func: x_{k+1} = f_func(x_k,t) where x_k is the state. The
%	function's second argument is time t for cases when the function
%	changes with time.
%
%	jacobian_func: the jacobian of f_func at state x and time t
%	jacobian_func(x,t)
%
%	t: current time
%
%	X_0: estimate of x at start_time
%
%	P_0_sqrt: square root factor of covariance P at
%	start_time. P_0 = P_0_sqrt*(P_0_sqrt')
%
%	Q_root: square root of process noise covariance matrix 
%	Q where Q = Q_root*(Q_root');
%
%OUTPUT:
%	estimate: estimate of x after predict phase 
%
%	covariance_sqrt: square root of covariance P after predict phase
%	P = covariance_sqrt*(covariance_sqrt')
	
	x = x_0;
	state_count = length(x);
	
	estimate = f_func(x,t);
	jac = jacobian_func(x,t);
	[~, covariance_sqrt] = qr([jac*P_0_sqrt,Q_root]');
	covariance_sqrt = covariance_sqrt';
	covariance_sqrt = covariance_sqrt(1:state_count,1:state_count);
end