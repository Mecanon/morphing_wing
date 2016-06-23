% Before running code, remember to check if csv files are in correct
% format, otherwise just open it in excel and save it.

x = zeros(15);

x(1) = 60E9;
x(2) = 60E9;
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
x(3) = 350.; %M_s
x(4) = 330;  %M_f
x(5) = 350.; %A_s
x(6) = 370.; %A_f

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
x(7) = 7.8E6;
x(8) = 7.3e6;

% Maximum and minimum transformation strain
x(9) = 0.09;
x(10) = 0.11;

x(11) = 0.021E-6;
% sigma_crit = 140E6;

% Coefficient of thermal expansion
% P.alpha = 0; %1E-5;

% Smoothn hardening parameters 
% NOTE: smoothness parameters must be 1 for explicit integration scheme
x(12) = 0.06; %0.618;
x(13) = 0.06; %0.313;
x(14) = 0.06; %0.759;
x(15) = 0.06; %0.358;

r = cost(x)