function output = cost(x)
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

%Read data from experiments. For constant stress sigma:
% - data_sigma(1): Temperature (in Celsius)
% - data_sigma(2): Strain
% - data_sigma(3): Stress

data_100 = textread('filtered_data_100MPa.txt');
T_100 = data_100(:,1) + 273.15;
eps_100 = data_100(:,2);
sigma_100 = data_100(:,3);

% INPUT:
% MATERIAL PARAMETERS (Structure: P)
% Young's Modulus for Austenite and Martensite 
P.E_A = x(1);
P.E_M = x(2);
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
P.M_s = x(3);
P.M_f = x(3) - x(4);
P.A_s = x(5);
P.A_f = x(6) + x(5);

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
P.C_A = x(7);
P.C_M = x(8);

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

for i = 1:(size(sigma_100,1))
    if sigma_100(i) > 100
        sigma_100(i:size(sigma_100,1)) = 100;
        break
    end
end
%Transform into MPa
sigma_100 = 1e6 * sigma_100;
% Elastic Prediction Check
elastic_check = 'N';

% Integration Scheme
integration_scheme = 'I';

%  try
[eps_num, MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T_100, sigma_100, P, elastic_check, integration_scheme );

figure()
box on
hold on
plot(T_100, eps_100,'b','LineWidth',1.5)
plot(T_100, eps_num,'r','LineWidth',1.5)

% Root-mean squared error:
output = sqrt(sum((eps_100-eps_num).^2)/numel(eps_100));
%  catch
% 	output = 1.;
end
