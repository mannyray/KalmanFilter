function [P] = sskf(A,C,Q,R);
%return steady covariance solution of continuous
%algebraic Riccati solution
%INPUT:
%	A: nxn matrix
%	C: mxn matrix (m <= n)
%	Q: nxn matrix
%	R: mxm matrix
%
%OUTPUT:
%	P: solution matrix P to continuous algebraic Riccati equation
	[P L G] = care(A',C',Q,R);
end
