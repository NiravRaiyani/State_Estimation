%% Function: mhe  (To perform constrained non-linear optimization for state estimation)
function X= mhe(Y,Xi,P)

global N ;

% Initial guess of states
X_init = [ones(N+1,1) ones(N+1,1)]';
X_init(:,1) = Xi;

% Defining the Uper and Lower bounds for the states
LB = [-10*ones(N+1,1)   -30*ones(N+1,1)]';
UB = [70*ones(N+1,1)  40*ones(N+1,1)]';

options = optimset('LargeScale','off','MaxFunEvals', 10e4,'MaxIter',10e4,'Display','notify'); 

% Constrained non-linear optimization for state estimation
X = fmincon('MHE_cost',X_init,[],[],[],[],LB,UB,[],options,Y,Xi,P);

end
