clc;
clear all;
close all;

%% DT State Space Model of the system

A = [0.9925 0; 0 0.985];    % 2 states : V - Volume and rho - Density
B = [1 1; -22.33 44.66];
C = [1 0; 0 1];
D = [0 0; 0 0];

% Disturbance Properties
Q = [1 0; 0 1];         % Process noice covariance
R = [0.1 0; 0 0.1];     % Measurement noice covariance
I = eye(2);

%% Initial Condition

X(:,1) = [0 , 0];

%% Input generation

t = [0:200]';
u1 = sin(t/2);  % Input Flow - F1
u2 = sin(t/2);  % Input Flow - F2
u = [u1 u2]';

%% Process Noice
n = length(t);
rng default
W1 = sqrt(Q(1,1))*randn(n,1);
W2 = sqrt(Q(2,2))*randn(n,1);
W = 0.1*[W1  W2]';

% Measurement Noice
V1 = sqrt(R(1,1))*randn(n,1);
V2 = sqrt(R(2,2))*randn(n,1);
V = [V1  V2]';

%% Real plant with process and measurement noice
for i = 1:1:200
    
X(:,i+1) = A*X(:,i) + B*u(:,i) + V(:,i);

Y(:,i) = C*X(:,i) + W(:,i);     

end
%% KALMAN FILTER

% Initial estimates of P
P = [0.1 0; 0 1];
% Initial estimates for states
X1(:,1) = [0.3 , 1];


for i = 1:1:200
    
K = C*P*(C*P*C + R);                        % Computing Kalman Gain

X1(:,i) = X1(:,i) + K*(Y(:,i) - C*X(:,i));  % Update Estimate

P = (I - K*C)*P;                            % Update error Covariance
    
X1(:,i+1) = A*X1(:,i) + B*u(:,i);           % Prediction of states for next time step

Y1(:,i) = C*X1(:,i);     
end

%% Real plant without disturbances
X2(:,1) = [0 , 0];
for i = 1:1:200
    
X2(:,i+1) = A*X2(:,i) + B*u(:,i);

Y2(:,i) = C*X2(:,i) ;     

end

%% Ploting the Results
n = 200;

% For Volume
figure()
plot(1:n,Y(1,1:n),1:n,Y1(1,1:n),'.',1:n,Y2(1,1:n),'--');
legend('Process with disturbances','Kalman Filter','Actual Process')
xlabel('Time(s)')
ylabel('Volume(m^3)')

% For Density
figure()
plot(1:n,Y(2,1:n),1:n,Y1(2,1:n),'.',1:n,Y2(2,1:n),'--');
legend('Process with disturbances','Kalman Filter','Actual Process')
xlabel('Time(s)')
ylabel('\rho(kg/m3)')









