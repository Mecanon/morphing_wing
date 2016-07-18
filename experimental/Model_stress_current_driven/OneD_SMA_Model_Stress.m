%--------------------------------------------------------------------------
%
% 1D SMA TAMU MODEL, EXPLICIT AND IMPLICIT IMPLEMENTATION 
%
% SUPLEMENTAL CODE FOR:
% Analysis and Design of Shape Memory Alloy Morphing Structures
% - Dimitris Lagoudas, Edwin Peraza Hernandez, Darren Hartl
%
% XXXXXXXXX ADD REFERENTIAL INFORMATION HERE XXXXXXXX
%
%
% APPLICATIONS:
% - 1D uniaxial loading of an SMA component
%
% CONDITIONS:
% - Initial fully austenitic state
%
%
% Authors: Cullen Nauck, Edwin Peraza Hernandez
% Department of Aerospace Engineering, Texas A&M University
% 3141 TAMU, College Station, TX 77843-3141
%
%--------------------------------------------------------------------------

% Initializing environment
clear all; close all; clc

%--------------------------------------------------------------------------
% INPUTS
%--------------------------------------------------------------------------

% Maximum Number of increments n per loading step
n = 300; 

% INPUT:
% Temperature and strain at the start and at the ends of each loading step
% Linear increments strain and temperature loading step assumed
T_inp = [30+273.15; 140+273.15; 30+273.15];

sigma_inp = [0; 0; 0];

% MATERIAL PARAMETERS (Structure: P)
% Young's Modulus for Austenite and Martensite 
P.E_A = 3.7427e+10;
P.E_M = 8.8888e+10;
% Transformation temperatures (M:Martensite, A:
% Austenite), (s:start,f:final)
P.M_s = 363.5013;
P.M_f = 297.9735;
P.A_s = 324.6427;
P.A_f = 385.0014;

% Slopes of transformation boundarings into austenite (C_A) and
% martensite (C_M) at Calibration Stress 
P.C_A = 7.1986e+06;
P.C_M = 7.9498e+06;

% Maximum and minimum transformation strain
P.H_min = 0.0387;
P.H_sat = 0.0550;

P.k = 4.6849e-09;
P.sig_crit = 0;

% Coefficient of thermal expansion
P.alpha = 0;

% Ambient T
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

% Generate strain and temperature states at each increment
% T: Temperature
for i = 1:(size(T_inp,1)-1)
    if i == 1
        T = linspace(T_inp(i), T_inp(i+1), n)';
    else     
        T = [T; linspace(T_inp(i), T_inp(i+1),n)'];
    end
end

% eps: Strain
for i = 1:(size(sigma_inp,1)-1)
    if i == 1
        sigma = linspace(sigma_inp(i), sigma_inp(i+1), n)';
    else     
        sigma = [sigma; linspace(sigma_inp(i), sigma_inp(i+1),n)'];
    end
end

elastic_check = 'N';

integration_scheme = 'I';

[eps,MVF,eps_t,E,MVF_r,eps_t_r ] = Full_Model_stress( T, sigma, P, elastic_check, integration_scheme );

figure()
box on 
plot(T, eps,'b','LineWidth',1.5)
xlabel('Temperature (K)')
ylabel('Strain')
title('One D SMA Models')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))

figure()
box on 
plot(T, MVF,'b','LineWidth',1.5)
xlabel('Temperature (K)')
ylabel('MVF')

figure()
box on 
plot(T, eps_t,'b','LineWidth',1.5)
xlabel('Temperature (K)')
ylabel('eps_t')
