function [estimate, covariance_sqrt] = cdekf_predict_phase(f_func,jacobian_func,h,dt_between_measurements,start_time,P_0_sqrt,x_0,Q_root)
%runs predict portion of update of the extended kalman filter 
%[estimate, covariance_sqrt] cdekf_predict_phase(f_func,jacobian_func,h,dt_between_measurements,start_time,P_0_sqrt,x_0,Q_root)
%INPUT:
%	f_func: encodes relationship \dot{x} = f_func(x,t) for state x
%	and time t.
%
%	jacobian_func: the jacobian of f_func at state x and time t
%
%	h: fixed timestep used by RK4 to cover dt_between_measurements 
%	to evolve estimate as well covariance matrix. h must be less
%	than dt_between_measurements
%
%	dt_between_measurements: distance between incoming measurements
%	 
%	start_time: start time over which solving RK4
%	(finish time is start_time + dt_between_measurements)
%
%	x_0: estimate of x at start_time
%
%	P_0_sqrt: square root factor of covariance P at
%	start_time.	P_0 = P_0_sqrt*(P_0_sqrt')
%
%	Q_root: square root of process noise covariance matrix 
%	Q where Q = Q_root*(Q_root');
%
%OUTPUT:
%	estimate: estimate of x at dt_between_measurements + start_time
%
%	covariance_sqrt: square root of covariance P after the 
%	update phase. P = covariance_sqrt*(covariance_sqrt')

	current_time = start_time;
	finish_time = start_time + dt_between_measurements;
	x = x_0;
	state_count = length(x);

	Phi = eye(state_count);
	Phi_sum_root = zeros(state_count,state_count);
		
	while(current_time < finish_time)
		if(current_time + h > finish_time)
			h = finish_time - current_time;
		end

		%RK4 for estimate
		k1 = h.*f_func(x,current_time);
		k2 = h.*f_func(x + k1./2,current_time+h/2);
		k3 = h.*f_func(x + k2./2,current_time+h/2);
		k4 = h.*f_func(x + k3, current_time + h);
		x = x + k1./6 + k2./3 + k3./3 + k4./6;
		
		%RK4 for phi
		jac_a  = jacobian_func(x,current_time);
		%k_phi_1 = h.*jac_A*(Phi);
		%k_phi_2 = k_phi_1 + (0.5*h).*jac_A*(k_phi_1);
		%k_phi_3 = k_phi_1 + (0.5*h).*jac_A*(k_phi_2);
		%k_phi_4 = k_phi_1 + h.*jac_A*(k_phi_3);
		%Phi_next = Phi + k_phi_1./6 + k_phi_2./3 + k_phi_3./3 + 
		%k_phi_4./6;
		%since jac_A is constant at current time the commented lines
		%can be replaced with the stuff below
		%with th
		%following few lines
		jac_a_2 = jac_a*jac_a;
		jac_a_3 = jac_a_2*jac_a;
		jac_a_4 = jac_a_3*jac_a;
		Phi_next = (eye(state_count) + (h).*jac_a + (1/2*h^2).*jac_a_2 +...
            (1/6*h^3).*jac_a_3 + (1/24*h^4).*jac_a_4)*Phi;

		%Trapezoid rule. add phi and phi_next together
		[~, R_tmp1] = qr([Phi*Q_root,Phi_next*Q_root]');
		R_tmp1 = sqrt(h*0.5).*R_tmp1';
		R_tmp1 = R_tmp1(1:state_count, 1:state_count);
		
		%Add previous sum to overall phi sum
		[~, R_tmp1] = qr([R_tmp1,Phi_sum_root]'); 
		R_tmp1 = R_tmp1';
		R_tmp1 = R_tmp1(1:state_count, 1:state_count);
		Phi_sum_root = R_tmp1;

		Phi = Phi_next;
		current_time = current_time + h;
	end
	
	[~, R_tmp1] = qr([Phi*P_0_sqrt,Phi_sum_root]');
	covariance_sqrt = R_tmp1';
	covariance_sqrt = covariance_sqrt(1:state_count,1:state_count);
	
	estimate = x;
end
