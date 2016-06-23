% Inputs:
% - x(1): E_M
% - x(2): E_A
% - x(3): M_s
% - x(4): M_s - M_f
% - x(5): A_s
% - x(6): A_f - A_s
% - x(7): C_M
% - x(8): C_A
% - x(9): H_min
% - x(10): H_max - H_min
% - x(11): k
% - x(12): n_1 
% - x(13): n_2
% - x(14): n_3
% - x(15): n_4
%alphas and sigma_crit are equal to zero in this problem
% Initial guess (parameters calibrated from experiments)
fun = @(x)cost(x);
% x0 = [37e9, 53e9, ...
%      273.15 + 55.74, 55.74 - 35.57, 273.15 + 72.93, 87.77 - 72.93, 4.54E6, 8.64E6, ...
%      0.1, 0.12, 0.0001, ...
%      .5, .5, .5, .5];
x0 = [60e9, 60e9, ...
     374.2813, 374.2813 - 290.7571, 318.9352,  397.9732 - 318.9352, 7.8E6, 7.8E6, ...
     0.0952, -0.001, 0.0001, ...
     .2215, .2059, .2040, .2856];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [1e9, 1e9, ...
     273.15, 0, 273.15, 0, 1E6, 1E6, ...
     0., -0.1, 0, ...
     0., 0., 0., 0.];
ub = [1e11, 1e11, ...
     398.15, 100., 398.15, 100., 10E6, 10E6, ...
     0.2, 0.3, 0.01, ...
     1., 1., 1., 1.];
nonlcon = [];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'MaxFunEvals', 100000);
x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

plot_optimized(x)

