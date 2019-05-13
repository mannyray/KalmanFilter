function dXdt = mRiccati(t, X, A, C,R, Q)
X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
dXdt = X*(A') + A*X - X*(C')*(inv(R))*C*X + Q; %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
%sourced from https://www.mathworks.com/matlabcentral/answers/305484-how-to-solve-a-riccati-control-differential-equation
