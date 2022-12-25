clear variables;
clc;

%Values of system variables
M = 1000;   % kg
m1 = 100;   % kg
m2 = 100;   % kg
l1 = 20;    % m
l2 = 10;    % m
g = 9.8;    % m/s^2

% System dynamics matrix
a = [0 1 0 0 0 0;
     0 0 -g*m1/M 0 -g*m2/M 0;
     0 0 0 1 0 0;
     0 0 -(M*g + m1*g)/(M*l1) 0 -m2*g/(M*l1) 0;
     0 0 0 0 0 1;
     0 0 -m1*g/(M*l2) 0 -(M*g + m2*g)/(M*l2) 0];
eigen_a = eig(a);