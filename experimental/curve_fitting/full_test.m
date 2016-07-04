% Before running code, remember to check if csv files are in correct
% format, otherwise just open it in excel and save it.
global initial_error
global initial_delta_eps
initial_error = 0;
initial_delta_eps = 0;
x = zeros(1,15);

x(1) = 38.03E9;
x(2) = 18.46E9;
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
x(3) = 273.15 + 62.38; %M_s
x(4) = 62.38 - 51.69;  %M_s - M_f
x(5) = 273.15 + 70.96; %A_s
x(6) = 83.49 - 70.96; %A_f - A_s

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
x(7) = 7.64E6; %C_M
x(8) = 8.15e6; %C_A

% Maximum and minimum transformation strain
x(9) =  0.0924617829295; %H_min
x(10) = 0.126325797371 -  0.0924617829295; %H_max - H_min

x(11) = 0.00524735758484e-6; %k
% sigma_crit = 140E6;

% Coefficient of thermal expansion
% P.alpha = 0; %1E-5;

% Smoothn hardening parameters 
% NOTE: smoothness parameters must be 1 for explicit integration scheme
x(12) = 0.6; %0.618;
x(13) = 0.6; %0.313;
x(14) = 0.6; %0.759;
x(15) = 0.6; %0.358;

r = full_cost(x)

% full_plot_optimized(x)