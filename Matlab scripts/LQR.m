clear variables;
clc;

syms L1 L2 m1 m2 g  M 

L1 = 20;    % m
L2 = 10;    % m
m1 = 100;   % kg
m2 = 100;   % kg
g = 9.8;    % m/s^2
M = 1000;   % kg

% System dynamics
A = [0 1 0 0 0 0;
     0 0 -g*m1/M 0 -g*m2/M 0;
     0 0 0 1 0 0;
     0 0 -(M*g + m1*g)/(M*l1) 0 -m2*g/(M*l1) 0;
     0 0 0 0 0 1;
     0 0 -m1*g/(M*l2) 0 -(M*g + m2*g)/(M*l2) 0];

B = transpose([0 1/M 0 1/(l1*M) 0 1/(l2*M)]);
% Output is described by variables x(t),theta1(t),theta2(t).
C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];
% Force is the only input parameter
D = transpose([1 0 0]);

% LQR weights
Q = [10 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 2500 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 2500 0;
     0 0 0 0 0 0];
R = 0.001;

% Get LQR parameters
[K, S, P] = lqr(A, B, Q, R);

% Calculating eigen values after applying LQR
% and getting K values, to check stability using
% Lyapunov's indirect method.
eigen_A = eig(A);
A_LQR = A-B*K;
eigen_A_LQR = eig(A_LQR);

% Define linear system model for LQR controller
system = ss(A-B*K, B, C, D);
tspan = 0:0.1:100;
initial = [0 0 deg2rad(10) 0 deg2rad(27) 0];
% Get system response for linear model 
[t,q1] = ode45(@(t,Q) linear_model(t,Q,-K*Q),tspan,initial);
figure(1);
hold on
plot(t,q1(:,1),'r')
plot(t,q1(:,3),'g')
plot(t,q1(:,5),'b')
ylabel('x(t),theta1(t),theta2(t)')
xlabel('time')
title('System response for Linearized system')
legend('x(t)','theta1(t)','theta2(t)')

% Get system response for non-linear model
[t,q2] = ode45(@(t,Q)nonLinear_model(t,Q,-K*Q),tspan,initial);
figure(2);
hold on
plot(t,q2(:,1),'r')
plot(t,q2(:,3),'g')
plot(t,q2(:,5),'b')
ylabel('x(t),theta1(t),theta2(t)')
xlabel('time')
title('System Response for the Non Linearized System')
legend('x(t)','theta1(t)','theta2(t)')