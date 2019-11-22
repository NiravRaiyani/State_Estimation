clc;
clear all;
close all;

%% Simulation and Sampling time 

T = 50; % Simulation time
ts= 0.05;   % Sampling time(Ts)

% 2 states : V - Volume and rho - Density
C = [1 0; 0 1];
D = [0 0; 0 0];

%% Disturbance Properties
Q = [0.05 0; 0 0.1];         % Process noice covariance
R = [10 0; 0 5];     % Measurement noice covariance
I = eye(2);
%% Initial Condition

X(:,1) = [20 , 845.333];

%% Input generation

t = [0:ts:T]';
u1 = 0.6*sin(t/2);  % Input Flow - F1
u2 = 0.5*sin(t/2);  % Input Flow - F2
u = [u1 u2]';

%% Process Noice
n = T/ts;
rng default
W1 = sqrt(Q(1,1))*randn(n,1);
W2 = sqrt(Q(2,2))*randn(n,1);
W = [0.5*W1  2*W2]';

% Measurement Noice
V1 = sqrt(R(1,1))*randn(n,1);
V2 = sqrt(R(2,2))*randn(n,1);
V = [V1  2*V2]';

%waitbar(0.1,h,'10% done')
%% Measurement data generation

for i = 1:1:(T/ts -1)
                                                                            %Addition of Process noice 'W'.
X(1,i+1) = X(1,i) +  ts*( -0.015*sqrt(X(1,i)) + u(1,i) + u(2,i)) + W(1,i);

X(2,i+1) = X(2,i) + (ts/X(1,i))*(((823 - X(2,i))*u(1,i)) + ((890 - X(2,i))*u(2,i))) + W(2,i) ;
                                                                            %Addition of Measurement noice 'V'.
Y(:,i+1) = C*X(:,i+1) + V(:,i+1);     
Y(:,1)   = C*X(:,1) + V(:,1);
                                                                            % Measurement with out any noice for comparison of results.
Y2(:,i+1) = C*X(:,i+1); 
Y2(:,1)   = C*X(:,1);
end


%% Extanded Kalman Filter

% Initial estimates of P
P = [0.5 0; 0 1.5];

% Initial estimates for states
X1(:,1) = [18 , 842];


for i = 1:1:(T/ts -1)
    
X1(1,i+1) = X1(1,i) +  ts*( -0.015*sqrt(X1(1,i)) + u(1,i+1) + u(2,i+1));  % Prediction of states for next time step                   
                                                    
X1(2,i+1) = X1(2,i) + (ts/X1(1,i))*(((823 - X1(2,i))*u(1,i+1)) + ((890 - X1(2,i))*u(2,i+1)));   
 
A = lin1(X1(:,i),u(:,i+1));                                           % Jacobian Evaluation at predicted states
    
P = A*P*A' + Q;                                                     % Covariance Prediction

N(:,i+1) = Y(:,i+1) - C*X1(:,i+1);                                  % Residual Measurement

S = C*P*C' + R;                                                     % Residual Covariance
    
K = C*P*inv(C*P*C' + R);                                            % Kalman Gain                            

X1(:,i+1) = X1(:,i+1) + K*N(:,i);                                   % Update estimate

Y1(:,i+1) = C*X1(:,i+1);                                            % Update Measurement

P = (I - K*C)*P;                                                    % Update Covariance
    
                                                                    
end
Y1(:,1)   = C*X1(:,1) + V(:,1);


%% Ploting the Results
n = T/ts ;

% For Volume
figure()
plot(0:ts:(T-ts),Y(1,1:n),'LineWidth',0.8);
hold on;
plot(0:ts:(T-ts),Y2(1,1:n),'--k',0:ts:(T-ts),Y1(1,1:n),'-r','LineWidth',1)
legend('Process with measurement noice','Process W/O measurement noice','Kalman Filter')
xlabel('Time(s)')
ylabel('Volume(m^3)')

% For Density
figure()
plot(0:ts:(T-ts),Y(2,1:n),'LineWidth',0.8);
hold on;
plot(0:ts:(T-ts),Y2(2,1:n),'--k',0:ts:(T-ts),Y1(2,1:n),'-r','LineWidth',1)
legend('Process with measurement noice','Process W/O measurement noice','Kalman Filter')
xlabel('Time(s)')
ylabel('\rho(kg/m3)')



%% Function : lin1 (to compute Jacobian at different states)

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









