clear all

%% Defining variables
syms m1 g m2 M L1 L2
m1 = 100;
m2 = 100;
M = 1000;
L1 = 20;
L2 = 10;
g = 9.81;
initial_states = [5 0 deg2rad(10) 0 deg2rad(35) 0];
tspan = 0:0.3:10;

%% Observability Check
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0;
    0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1;
    0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
C2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = [0; 0; 0];
Observability1 = rank([C1' A'*C1' ((A')^2)*C1' ((A')^3)*C1' ((A')^4)*C1' ...
    ((A')^5)*C1']);
Observability2 = rank([C2' A'*C2' ((A')^2)*C2' ((A')^3)*C2' ((A')^4)*C2' ...
    ((A')^5)*C2']);
Observability3 = rank([C3' A'*C3' ((A')^2)*C3' ((A')^3)*C3' ((A')^4)*C3' ...
    ((A')^5)*C3']);
Observability4 = rank([C4' A'*C4' ((A')^2)*C4' ((A')^3)*C4' ((A')^4)*C4' ...
    ((A')^5)*C4']);
% We note that, for state vector ((theta1(t),theta2(t)), the system is not
% observable. Hence, we wouldn't design an observer for that.

Actual_sys1 = ss(A,B,C1,D);
Actual_sys3 = ss(A,B,C3,D);
Actual_sys4 = ss(A,B,C4,D);


%% Luenberger Observer Design
Process_noise = 0.5*eye(6);                %Process Noise
Measurement_noise = 0.1*eye(3);           %Measurement Noise
[L_Matrix1,~,~] = lqe(A,Process_noise,C1,Process_noise,Measurement_noise);
[L_Matrix3,~,~] = lqe(A,Process_noise,C3,Process_noise,Measurement_noise);
[L_Matrix4,P,E] = lqe(A,Process_noise,C4,Process_noise,Measurement_noise);

A_Closedloop1 = A-(L_Matrix1*C1);
A_CLosedloop3 = A-(L_Matrix3*C3);
A_CLosedloop4 = A-(L_Matrix4*C4);

Estimated_sys1 = ss(A_Closedloop1,[B L_Matrix1],C1,0);
Estimated_sys3 = ss(A_CLosedloop3,[B L_Matrix3],C3,0);
Estimated_sys4 = ss(A_CLosedloop4,[B L_Matrix4],C4,0);

%% Estimated vs actual state vectors for step input
unitStep(10:length(tspan)) = 0.5;

[y1,~] = lsim(Actual_sys1,unitStep,tspan);
[x1,~] = lsim(Estimated_sys1,[unitStep;y1'],tspan);

[y3,~] = lsim(Actual_sys3,unitStep,tspan);
[x3,~] = lsim(Estimated_sys3,[unitStep;y3'],tspan);

[y4,~] = lsim(Actual_sys4,unitStep,tspan);
[x4,t] = lsim(Estimated_sys4,[unitStep;y4'],tspan);

figure();
hold on
plot(t,y1(:,1),'r','Linewidth',2)
plot(t,x1(:,1),'k--','Linewidth',1)
ylabel('States')
xlabel('time(sec)')
legend('Actual x(t)','Estimated x(t)')
title('Step input and Response for output vector: (x(t))')
hold off

figure();
hold on
plot(t,y3(:,1),'k','Linewidth',2)
plot(t,y3(:,3),'m','Linewidth',2)
plot(t,x3(:,1),'r--','Linewidth',1)
plot(t,x3(:,3),'b--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Step input and Response for output vector: (x(t),theta2(t))')
hold off

figure();
hold on
plot(t,y4(:,1),'k','Linewidth',2)
plot(t,y4(:,2),'r','Linewidth',2)
plot(t,y4(:,3),'b','Linewidth',2)
plot(t,x4(:,1),'m--','Linewidth',1)
plot(t,x4(:,2),'g--','Linewidth',1)
plot(t,x4(:,3),'b--','Linewidth',1)
ylabel('State Variables')
xlabel('t')
legend('x(t)','theta1(t)','theta2(t)','Estimated x(t)', ...
    'Estimated theta1(t)','Estimated theta2(t)')
title(['Step input and Response for output vector:' ...
    ' (x(t),theta1(t),theta2(t))'])
hold off

%% Observer design response for Linearized system
[t,q1] = ode45(@(t,q)linear_model_Observer1(t,q,L_Matrix1), ...
    tspan,initial_states);
figure();
hold on
plot(t,q1(:,1))
ylabel('states')
xlabel('t)')
title('Linear Observer for output vector: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)linear_model_Observer3(t,q,L_Matrix3),tspan, ...
    initial_states);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('states')
xlabel('t')
title('Linear Observer for output vector: (x(t),theta_2(t))')
legend('x','theta2')
hold off

[t,q4] = ode45(@(t,q)linear_model_Observer4(t,q,L_Matrix4),tspan, ...
    initial_states);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('states')
xlabel('t')
title(' Linear Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

%%
%% Observer design response for Non Linear system
[t,q1] = ode45(@(t,q)nonLinear_model_Observer1(t,q,1,L_Matrix1),tspan, ...
    initial_states);
figure();
hold on
plot(t,q1(:,1))
ylabel('states')
xlabel('t)')
title('Non-Linear System Observer for output vector: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)nonLinear_model_Observer3(t,q,1,L_Matrix3),tspan, ...
    initial_states);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('states')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)nonLinear_model_Observer4(t,q,1,L_Matrix4),tspan, ...
    initial_states);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('states')
xlabel('t')
title(['Non-Linear System Observer for output vector: (x(t),theta_1(t),' ...
    'theta_2(t))'])
legend('x','theta_1','theta_2')
hold off