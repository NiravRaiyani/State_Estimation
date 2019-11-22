clc;
clear all;
close all;
f = waitbar(0,'Please wait...');
pause(.5)

global N C ts Q R A B;
T  = 90;                                            % Simulation Time
ts = 1;                                             % Sampling Time
N  = 10;                                            % Moving Horizon Window length

A = [0.99 0.2; -0.1 0.3];
B = [1 1]';


% 2 measurements : H - Level in the tank and rho - Density
C = [1 -3];
D = 1;

% Initial estimates of P
P = [1 0; 0 1];

%% Property of process and measurement noise
Q = 1;                                     % Process Noice co - variance
R = 0.01;                                  % Measurement Noice Co - variance

%% Input Generation
t = [0:ts:T]';
u1 = 0.6*sin(t/2);                         % Input Flow - F1
u2 = 0.5*sin(t/2);                         % Input Flow - F2
u = [u1 u2]'; 

%% Noise generation
n = T/ts;
rng default
% Process Noice
W = sqrt(Q(1,1))*randn(n,1);
W = W';

% Measurement Noice
V = sqrt(R(1,1))*randn(n,1);
V = V';
%% Process Output Measurement

% Initial Condition
X(:,1) = [0 , 0];

for i = 1:1:(T/ts)

X(:,i+1) =  A*X(:,i)+ B*W(:,i)
    
    
Y(:,i) = C*X(:,i)+ V(:,i);     

end

%% MHE - Moving Horizon Estimation
Xm(:,1) = [0 , 0];
for i = N+1:1:T/ts           
     
    
   
    Xm(:,i-N:i) = mhe(Y(:,i-N:i-1),Xm(:,i-N),P(:,:,i-N));
    Ym(:,i-N:i) = C*Xm(:,i-N:i);
    
    
    %To calculate the Jacobian for given X and U
%     A(:,:,i-N+1) = lin1(Xm(:,i-N+1),u(:,i-N+1));
    
    % The Riccati equation to calculate P(k+1)
    P(:,:,i-N+1) = A*P(:,:,i-N)*A' - A*P(:,:,i-N)*C'*inv(C*P(:,:,i-N)*C' + R)*C*P(:,:,i-N)*A' + Q;

    
    
    perc = round((i/(T/ts))*100);                                   
    waitbar(i/(T/ts),f,sprintf('Calculating states by MHE: %d%% done!!!',perc));
end

%% Plots of the results
n = T/ts ;
%subplot(2,1,1)
% For Volume
figure()
plot(0:ts:(T-ts),X(1,1:n),'LineWidth',0.8);
hold on;
plot(0:ts:(T-ts),Xm(1,1:n),'-r','LineWidth',1.5)
legend('True','MHE')
xlabel('Time(s)')
ylabel('X_1')

%subplot(2,1,2)
% For Density
figure()
plot(0:ts:(T-ts),X(2,1:n),'LineWidth',0.8);
hold on;
plot(0:ts:(T-ts),Xm(2,1:n),'-r','LineWidth',1.5)
legend('True','MHE')
xlabel('Time(s)')
ylabel('X_2')


