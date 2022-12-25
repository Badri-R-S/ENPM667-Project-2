clear all

%% Defining variables
syms m1 g m2 M L1 L2 x dx 
m1 = 100;
m2 = 100;
M = 1000;
L1 = 20;
L2 = 10;
g = 9.81;
tspan = 0:0.1:100;
% q = [x dx t1 dt1 t2 dt2];
%Enter initial conditions
initial_states = [2 0 deg2rad(2) 0 deg2rad(5) 0];

%% Linearized Model
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0;
    0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1;
    0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
D = [1;0;0];
Actual_sys1 = ss(A,B,C1,D);

%% LQR Controller
Q = [100 0 0 0 0 0;
     0 1000 0 0 0 0;
     0 0 3000 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 3000 0;
     0 0 0 0 0 3000];
R = 0.001;
[K,S,~] = lqr(A,B,Q,R);
LQR_sys = ss(A-B*K,B,C1,D);
% step(sys,200);

%% Kalman Estimator Design
Process_noise = 0.01*eye(6);                %Process Noise
Measurement_noise = 0.001;                      %Measurement Noise
[L_Matrix1,P,E] = lqe(A,Process_noise,C1,Process_noise ...
    ,Measurement_noise*eye(3)); %Considering vector output: x(t)
Ac1 = A-(L_Matrix1*C1);
Estimated_sys1 = ss(Ac1,[B L_Matrix1],C1,0);

%% Non-linear Model LQG Response
[t,q1] = ode45(@(t,q)nonLinear_model_Observer1(t,q,-K*q,L_Matrix1), ...
    tspan,initial_states);
figure();
hold on
plot(t,q1(:,1))
ylabel('state')
xlabel('time')
title('LQG output for the Non-Linear system with state vector: x(t)')
legend('x')
hold off