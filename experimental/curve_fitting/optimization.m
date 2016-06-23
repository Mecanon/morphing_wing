% Inputs:
% - x(1): E_M
% - x(2): E_A
% - x(3): M_s
% - x(4): M_f
% - x(5): A_s
% - x(6): A_f
% - x(7): C_M
% - x(8): C_A
% - x(9): H_min
% - x(10): H_sat
% - x(11): k
% - x(12): n_1 
% - x(13): n_2
% - x(14): n_3
% - x(15): n_4
%alphas and sigma_crit are equal to zero in this problem
% Initial guess (parameters calibrated from experiments)
fun = @(x)cost(x);
x0 = [37e9, 53e9, ...
     350.0, 330.0, 350.0, 370.0, 7.8E6, 7.8E6, ...
     0.1, 0.12, 0.0001, ...
     .5, .5, .5, .5];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [1e9, 1e9, ...
     273.15, 273.15, 273.15, 273.15, 1E6, 1E6, ...
     0., 0.05, 0., ...
     0., 0., 0., 0.];
ub = [1e11, 1e11, ...
     398.15, 398.15, 398.15, 398.15, 10E6, 10E6, ...
     0.2, 0.2, 0.01, ...
     1., 1., 1., 1.];
nonlcon = [];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'MaxFunEvals', 10000);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);