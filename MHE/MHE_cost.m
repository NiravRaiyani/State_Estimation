%% Function: MHE_cost (To calculate the total cost at every time step)
function J = MHE_cost(X,Y,Xi,P)

global N C Q R;

% Storing the values of the current states
Xc   = X;
Yc   = Y;

X_ic =Xi;

% Initial cost
J = 0;

% Arrival cost
IC = (X_ic - Xc(:,1))'*inv(P)*(X_ic - Xc(:,1));

% Cost function Evaluation
for k = 1:1:N
    
V(:,k) = Yc(:,k) - C*Xc(:,k);
W(:,k) = Xc(:,k+1) - sys(Xc(:,k));

J = J + V(:,k)'*inv(R)*V(:,k) + W(:,k)'*inv(Q)*W(:,k);

end


J = J + IC;

end