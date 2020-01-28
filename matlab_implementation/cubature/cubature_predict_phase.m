function [estimate, covariance_sqrt] = cubature_predict_phase(f_func,t,P_0_sqrt,x_0,Q_root)
%runs predict portion of update of the dd-ekf
%[estimate, covariance_sqrt] cubature_predict_phase(f_func,t,P_0_sqrt,x_0,Q_root)
%INPUT:
%	f_func: x_{k+1} = f_func(x_k,t) where x_k is the state. The
%	function's second argument is time t for cases when the function
%	changes with time.
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
	
	sigma_points = zeros(state_count,state_count*2);
	sigma_points(:,1:state_count) = x_0 + P_0_sqrt*sqrt(state_count).*eye(state_count,state_count);
	sigma_points(:,(state_count+1):end) = x_0 - P_0_sqrt*sqrt(state_count).*eye(state_count,state_count);
	
	for i=1:(2*state_count)
		sigma_points(:,i) = f_func(sigma_points(:,i),t);
	end
	
	estimate = (1/(2*state_count))*sum(sigma_points,2);
	
	covariance = zeros(state_count,state_count);
	for i=1:(2*state_count)
		covariance = covariance + (sigma_points(:,i)-estimate)*(sigma_points(:,i)-estimate)';
	end
	covariance = 1/(2*state_count) * covariance + Q_root*Q_root';
	
	covariance_sqrt = chol(covariance)';
end
