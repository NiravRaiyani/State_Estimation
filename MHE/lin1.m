function P1 = lin1(X,u)

% Predicted Estimates
Vs = X(1,1); % Volume Estimate
rs = X(2,1); % Density Estimate
U1 = u(1,1); % Input : Flow F1
U2 = u(2,1); % Input : Flow F2



% Jacobian 
A = [(-0.0075*(1/sqrt(Vs)))  0; (-1)*(1/(Vs^2))*(((823-rs)*U1)+((890 - rs)*U2))  ((-U1/Vs) + (-U2/Vs))];
B = [1 1 ; ((823 - rs)/Vs)  ((890 - rs)/Vs)];
C = [1 0 ; 0 1];
D = [0 0 ; 0 0];

% CT Linearised States Space Model
sys_c = ss(A,B,C,D);

% DT Linearised State Space Model
sys_D = c2d(sys_c,0.1); % Discritisation with Ts = 1.

P1 = sys_D.A; % Returning the value of A matrix for covariance calculation
end


