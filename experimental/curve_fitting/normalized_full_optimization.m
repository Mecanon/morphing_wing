global initial_error
global initial_delta_eps
initial_error = 0;
initial_delta_eps = 0;

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

% x0 = [37e9, 53e9, ...
%      273.15 + 55.74, 55.74 - 35.57, 273.15 + 72.93, 87.77 - 72.93, 4.54E6, 8.64E6, ...
%      0.1, 0.12, 0.0001, ...
%      .5, .5, .5, .5];
% x0 = [60e9, 60e9, ...
%      374.2813, 374.2813 - 290.7571, 318.9352,  397.9732 - 318.9352, 7.8E6, 7.8E6, ...
%      0.0952, -0.001, 0.0001, ...
%      .2215, .2059, .2040, .2856];
x0 = [38.03E9, 38.03E9 - 18.46E9, 273.15 + 62.38, 62.38 - 51.69, ...
     273.15 + 70.96, 83.49 - 70.96, 8.15e6, 7.64E6,...
     0.0937, 0.1275 - 0.0937, 0.00524e-6, ...
     0.2, 0.2, 0.2, 0.3];

A = [];
b = [];
Aeq = [];
beq = [];
lb = [20e9, 5e9, ...
     273.15 + 50, 0, 273.15 + 40, 0, 4E6, 4E6, ...
     0.05, 0., 0.001e-6, ...
     0., 0., 0., 0.];
 
ub = [60e9, 40e9, ...
     273.15 + 140, 70., 273.15 + 140, 70., 10E6, 10E6, ...
     0.15, 0.05, 0.01e-6, ...
     1., 1., 1., 1.];

% Normalized x0
n_x0 = (x0 - lb)./(ub-lb);

% Normalized lower and upper bounds
n_lb = zeros(size(lb));
n_ub = ones(size(ub));

% Define function to be optimized
fun = @(x)normalized_full_cost(x, lb, ub);

nonlcon = [];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'MaxFunEvals', 100000, 'PlotFcns',{@optimplotx,...
    @optimplotfval,@optimplotfirstorderopt});
x = fmincon(fun, n_x0, A, b, Aeq, beq, n_lb, n_ub, nonlcon, options);
% opts = gaoptimset(...
%     'PopulationSize', 20, ...
%     'Generations', 50, ...
%     'Display', 'iter', ...
%     'EliteCount', 2);
% x = ga(fun, 15, A, b, Aeq, beq, lb, ub, nonlcon, opts);

normalized_full_plot_optimized(x, lb, ub)

