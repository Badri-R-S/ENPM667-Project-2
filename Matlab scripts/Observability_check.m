%% Defining variables
syms L1 L2 m1 m2 g  M 

% Defining system dynamics
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0;
    0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1;
    0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C_matrix_1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
C_matrix_2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
C_matrix_3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
C_matrix_4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
%% Observability Check
Observability_1 = rank([C_matrix_1' A'*C_matrix_1' ((A')^2)*C_matrix_1' ...
    ((A')^3)*C_matrix_1' ((A')^4)*C_matrix_1' ((A')^5)*C_matrix_1']);
Observability_2 = rank([C_matrix_2' A'*C_matrix_2' ((A')^2)*C_matrix_2' ...
    ((A')^3)*C_matrix_2' ((A')^4)*C_matrix_2' ((A')^5)*C_matrix_2']);
Observability_3 = rank([C_matrix_3' A'*C_matrix_3' ((A')^2)*C_matrix_3' ...
    ((A')^3)*C_matrix_3' ((A')^4)*C_matrix_3' ((A')^5)*C_matrix_3']);
Observability_4 = rank([C_matrix_4' A'*C_matrix_4' ((A')^2)*C_matrix_4' ...
    ((A')^3)*C_matrix_4' ((A')^4)*C_matrix_4' ((A')^5)*C_matrix_4']);