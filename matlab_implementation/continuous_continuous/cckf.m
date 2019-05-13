function [covariances] = cckf(t_start,t_final,outputs,P_0,A,C,R,Q)
%Solves Differential Riccati equation by using very simple 
%Euler scheme. Hence this implementation should not be used
%for higher mode count or anything 'serious'
%INPUT:
%	t_start:start time of solving Riccati equation
%
%	t_final: finish time
%
%	outputs: The amount of timesteps to be taken to solve
%	P(t)
%
%	P_0: P(t_start)
%
%	A: linear system matrix
%
%	C: observation matrix
%
%	R: sensor noise covariance matrix
%
%	Q: process noise covariance matrix
%OUTPUT:
%	covariances: cell of size outputsX1
%	where the ith entry of cell array is 
%	covariances{i,1} = P(t_start + (i-1)*dt)
%	where dt = (t_final-t_start)/outputs

dt = (t_final - t_start)/outputs;
covariances = cell(outputs,1);
current_time = t_start;

P = P_0;
counter = 1;
while current_time < t_final
	covariances{counter,1}  = P;

	k1 = reshape(dt.*mRiccati(current_time,P,A,C,R,Q),size(A));
	k2 = reshape(dt.*mRiccati(current_time + dt/2,P + 0.5.*k1,A,C,R,Q),size(A));
	k3 = reshape(dt.*mRiccati(current_time + dt/2,P + 0.5.*k2,A,C,R,Q),size(A));
	k4 = reshape(dt.*mRiccati(current_time + dt,P + k3,A,C,R,Q),size(A));

	P = (1/6).*(k1 + 2.*k2 + 2.*k3 + k4) + P;
	P = 0.5*(P + P');
	
	current_time = current_time + dt;
	counter = counter + 1;
end

end
