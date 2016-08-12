function output = normalized_full_plot_optimized(x, lb, ub)
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

% Denormalizing
x = x.*(ub - lb) + lb;

%Read data from experiments. For constant stress sigma:
% - data_sigma(1): Temperature (in Celsius)
% - data_sigma(2): Strain
% - data_sigma(3): Stress
data_50 = textread('filtered_data_50MPa.txt');
data_100 = textread('filtered_data_100MPa.txt');
data_150 = textread('filtered_data_150MPa.txt');
data_172 = textread('filtered_data_172MPa.txt');
data_200 = textread('filtered_data_200MPa.txt');

T_50 = data_50(:,1) + 273.15;
T_100 = data_100(:,1) + 273.15;
T_150 = data_150(:,1) + 273.15;
T_172 = data_172(:,1) + 273.15;
T_200 = data_100(:,1) + 273.15;

eps_50 = data_50(:,2);
eps_100 = data_100(:,2);
eps_150 = data_150(:,2);
eps_172 = data_172(:,2);
eps_200 = data_200(:,2);

sigma_50 = data_50(:,3);
sigma_100 = data_100(:,3);
sigma_150 = data_150(:,3);
sigma_172 = data_172(:,3);
sigma_200 = data_200(:,3);

% INPUT:
% MATERIAL PARAMETERS (Structure: P)
% Young's Modulus for Austenite and Martensite 
P.E_M = x(1);
P.E_A = x(1) - x(2);
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
P.M_s = x(3);
P.M_f = x(3) - x(4);
P.A_s = x(5);
P.A_f = x(6) + x(5);

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
P.C_M = x(7);
P.C_A = x(8);

% Maximum and minimum transformation strain
P.H_min = x(9);
P.H_sat = x(9) + x(10);

P.k = x(11);
P.sig_crit = 0;

% Coefficient of thermal expansion
P.alpha = 0.;

% Smoothn hardening parameters 
% NOTE: smoothness parameters must be 1 for explicit integration scheme
P.n1 = x(12);
P.n2 = x(13);
P.n3 = x(14);
P.n4 = x(15);

% Algorithmic delta for modified smooth hardening function
P.delta=1e-5;

% Calibration Stress
P.sig_cal=200E6;

% Tolerance for change in MVF during implicit iteration
P.MVF_tolerance=1e-8;

disp(P)
%Transform into MPa
sigma_50 = 1e6 * 50*ones(size(sigma_50));
sigma_100 = 1e6 * 100*ones(size(sigma_100));
sigma_150 = 1e6 * 150*ones(size(sigma_150));
sigma_172 = 1e6 * 172*ones(size(sigma_172));
sigma_200 = 1e6 * 200*ones(size(sigma_200));

% Elastic Prediction Check
elastic_check = 'N';

% Integration Scheme
integration_scheme = 'I';


[eps_num_50, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_50, sigma_50, P, elastic_check, integration_scheme );
[eps_num_100, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_100, sigma_100, P, elastic_check, integration_scheme );
[eps_num_150, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_150, sigma_150, P, elastic_check, integration_scheme );
[eps_num_172, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_172, sigma_172, P, elastic_check, integration_scheme );
[eps_num_200, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_200, sigma_200, P, elastic_check, integration_scheme );

eps_num_50 = eps_num_50 - min(eps_num_50);
eps_num_100 = eps_num_100 - min(eps_num_100);
eps_num_150 = eps_num_150 - min(eps_num_150);
eps_num_172 = eps_num_172 - min(eps_num_172);
eps_num_200 = eps_num_200 - min(eps_num_200);


figure()
box on
hold on
plot(T_50, eps_50,'b','LineWidth',1.5)
plot(T_50, eps_num_50,'r','LineWidth',1.5)
title('50MPa')

figure()
box on
hold on
plot(T_100, eps_100,'b','LineWidth',1.5)
plot(T_100, eps_num_100,'r','LineWidth',1.5)
title('100MPa')

figure()
box on
hold on
plot(T_150, eps_150,'b','LineWidth',1.5)
plot(T_150, eps_num_150,'r','LineWidth',1.5)
title('150MPa')

figure()
box on
hold on
plot(T_172, eps_172,'b','LineWidth',1.5)
plot(T_172, eps_num_172,'r','LineWidth',1.5)
title('172MPa')

figure()
box on
hold on
plot(T_200, eps_200,'b','LineWidth',1.5)
plot(T_200, eps_num_200,'r','LineWidth',1.5)
title('200MPa')
end