syms L1 L2 m1 m2 g  M 

L1 = 20;    % m
L2 = 10;    % m
m1 = 100;   % kg
m2 = 100;   % kg
g = 9.8;    % m/s^2
M = 1000;   % kg

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0;
    0 0 -((M*g)+(m1*g))/(M*L1)
    0 -g*m2/(M*L1) 0; 0 0 0 0 0 1; 
    0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
disp(C)
Rank = rank([B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B]);
Det = det([B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B]);