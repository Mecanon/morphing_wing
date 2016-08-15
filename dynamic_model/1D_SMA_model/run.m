
% Initializing environment
clear all; close all; clc

%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------

% Maximum Number of increments n per loading step
n = 100; 

% INPUT:
% Temperature and time at the start and at the ends of each loading step
% Linear increments strain and temperature loading step assumed
% data_ = textread('filtered_data_50MPa.txt');
% 
% t_inp = data(:,1);
% 
% T_inp = data(:,2) + 273.15;
% 
% eps_inp = data(:,3);

% Linear increments strain, temperature and current, using linspace function.
t = linspace(0,120,n)';
%eps_inp = linspace(0.10,0.04,4)';

% Pre-stress
sigma_o = 100e6; %Pa

% MATERIAL PARAMETERS (Structure: P)
% Young's Modulus for Austenite and Martensite 
P.E_A = 2.1496e+10;
P.E_M = 3.3453e+10;
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
P.M_s = 362.5851;
P.M_f = 297.4771;
P.A_s = 318.3625;
P.A_f = 386.8458;

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
P.C_A = 8036800;
P.C_M = 7123000;

% Maximum and minimum transformation strain
P.H_min = 0.0924;
P.H_sat = 0.1209;

P.k = 5.9713e-09;
P.sig_crit = 0;

%Energy Coefficients
% Coefficient of thermal expansion
P.alpha = 10E-6; % can be zero
% Mass density
P.rho= 6500; %kg/m^3
% Specific Heat
P.c= 837.36;
% Heat convection coefficient
P.h = 0; % 1 is True and 0 is False
% Wire electric resistivity
P.rho_E_M = 6e-6;
P.rho_E_A = 6e-6;
% Ambient Temperature (Initial Temperature??)
T_ambient=303.15; %K
% Model Geometry
% d: Diameter of considered 1D model
P.d = 0.381e-3;

% Smoothn hardening parameters 
P.n1 = 0.1919; %0.618;
P.n2 = 0.1823; %0.313;
P.n3 = 0.1623; %0.759;
P.n4 = 0.2188; %0.358;

% Algorithmic delta for modified smooth hardening function
P.delta=1e-5;

% Calibration Stress
P.sig_cal=200E6;

% Tolerance for change in MVF during implicit iteration
P.MVF_tolerance=1e-8;

% Generate strain and time states at each increment

% eps: Strain
eps = zeros((size(t,1)),1);
eps(1,1) = sigma_o/P.E_M + H_cursolver(sigma_o, P.sig_crit,P.k,P.H_min,P.H_sat);
% I: current
%I = zeros((size(t,1)),1);
%I(1,1) = 0;
I = linspace(0.,0.2,n)';
% Elastic Prediction Check
% prompt = {'Will the Elastic Prediction Check be Transformation Surface Rate-informed or not (Y/N)?'};
% dlg_title = '1D SMA Model Elastic Prediction Check';
% num_lines = 1;
% defaultans = {'N','hsv'};
% elastic_check = inputdlg(prompt,dlg_title,num_lines,defaultans);
elastic_check = 'N';
%%%%% Put the loop here

tol = 1e-6;
for i = 2:(size(t,1))
    error = 9999;
    counter = 0;
    while error>tol
        old_eps = eps(i-1,1);
        if i == 2
            eps(i,1) = sigma_o/P.E_M + H_cursolver(sigma_o, P.sig_crit,P.k,P.H_min,P.H_sat);
        else
            eps(i,1) = sigma(i-1,1)/E(i-1,1) + P.alpha*(T(i-1,1) - T_ambient) + eps_t(i-1,1);
            disp('sigma/E')
            disp(sigma(i-1,1)/E(i-1,1))
            disp('thermal')
            disp(P.alpha*(T(i-1,1) - T_ambient))
            disp('transf')
            disp(eps_t(i-1,1))
            disp('eps')
            disp(sigma(i-1,1)/E(i-1,1) + P.alpha*(T(i-1,1) - T_ambient) + eps_t(i-1,1))
        end
        [sigma,MVF,T,eps_t,E,MVF_r,eps_t_r ] = Full_Model_TC(i, t, eps, I, P, elastic_check, T_ambient, sigma_o);
        error = abs(eps(i,1) - old_eps);
        counter = counter + 1
    end
end
%[sigma,MVF,T,eps_t,E,MVF_r,eps_t_r ] = Full_Model_TC( t, eps, I, P, elastic_check, T_ambient );

figure(1)
box on 
plot(eps,sigma/1E6,'b','LineWidth',1.5)
xlabel('Strain')
ylabel('Stress (MPa)')
title('One D SMA Models')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))

figure(2)
box on 
plot(T,eps,'b','LineWidth',1.5)
xlabel('Temperature (K)')
ylabel('Strain')
title('One D SMA Models T vs. \epsilon')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))

figure(3)
box on 
plot(T,sigma/1E6,'b','LineWidth',1.5)
xlabel('Temperature (K)')
ylabel('Stress (MPa)')
title('One D SMA Models T vs. \sigma')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))

figure(4)
box on
hold on
plotyy(t,T,t,MVF);
ylabel('Temperature (K)')
xlabel('Time (s)')

figure(5)
box on
plot(t,sigma/1E6,'r','LineWidth',1.5)
ylabel('Stress (MPa)')
xlabel('Time (s)')

figure(6)
box on
plot(t,I,'r','LineWidth',1.5)
ylabel('Current')
xlabel('Time (s)')

figure(7)
box on
plot(t,eps_t,'r','LineWidth',1.5)
ylabel('Transformation strain')
xlabel('Time (s)')