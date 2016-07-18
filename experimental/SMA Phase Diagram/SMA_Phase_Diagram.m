%--------------------------------------------------------------------------
%
% SMA Stress-Temperature Phase Diagram
%
%
% SUPLEMENTAL CODE FOR:
% Analysis and Design of Shape Memory Alloy Morphing Structures
% - Dimitris Lagoudas, Edwin Peraza Hernandez, Darren Hartl
%
%
% XXXXXXXXX ADD REFERENTIAL INFORMATION HERE XXXXXXXX
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
%% INPUTS
%--------------------------------------------------------------------------

% Number of stress values
n = 200; 

% INPUT: Stress (sigma)

sigma_inp = [0; 200e6]; 

%% MATERIAL PARAMETERS (Structure: P)
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
P.sigma_crit = 0;

% Smoothn hardening parameters 
% NOTE: smoothness parameters must be 1 for explicit integration scheme
P.n1 = 0.1752; %0.618;
P.n2 = 0.1789; %0.313;
P.n3 = 0.1497; %0.759;
P.n4 = 0.2935; %0.358;

% Calibration Stress
P.sig_cal=200E6;

% Current transformation strain at calibration stress
H_cur_cal = H_cursolver(P.sig_cal, P.sigma_crit,P.k,P.H_min,P.H_sat);

% Partial Derivative of H_cur at calibration stress (dH_cur)
dH_cur=partial_Hcur_sigma(P.sig_cal,P.sigma_crit,P.k,P.H_sat,P.H_min);

%% Transformation Parameters (structure: TP)
TP.rho_delta_s0 = (-2*(P.C_M*P.C_A)*(H_cur_cal+P.sig_cal*dH_cur+P.sig_cal*(1/P.E_M-1/P.E_A)))/(P.C_M+P.C_A);
TP.D = ((P.C_M-P.C_A)*(H_cur_cal+P.sig_cal*dH_cur+P.sig_cal*(1/P.E_M-1/P.E_A)))/((P.C_M+P.C_A)*(H_cur_cal+...
    P.sig_cal*dH_cur));
TP.a1 = TP.rho_delta_s0*(P.M_f-P.M_s);
TP.a2 = TP.rho_delta_s0*(P.A_s-P.A_f);
TP.a3 = -TP.a1/4*(1+1/(P.n1+1)-1/(P.n2+1))+TP.a2/4*(1+1/(P.n3+1)-1/(P.n4+1));
TP.rho_delta_u0 = TP.rho_delta_s0/2*(P.M_s+P.A_f);
TP.Y_0_t = TP.rho_delta_s0/2*(P.M_s-P.A_f)-TP.a3;


%% Determine Stress vs. Temperature PHase Diagram

%Generate Stress (sigma) at each increment
for i = 1:(size(sigma_inp,1)-1)
    if i == 1
        sigma = linspace(sigma_inp(i), sigma_inp(i+1), n)';
    else     
        sigma = [sigma; linspace(sigma_inp(i), sigma_inp(i+1),n)'];
    end
end

% Arrays of output variables
% T_fwd_0: Temperature array for forward transformation at MVF=0
T_fwd_0 = zeros((size(sigma,1)),1);

% T_fwd_1: Temperature array for forward transformation at MVF=1
T_fwd_1 = zeros((size(sigma,1)),1);

% T_rev_0: Temperature array for reverse transformation at MVF=0
T_rev_0 = zeros((size(sigma,1)),1);

% T_rev_0: Temperature array for reverse transformation at MVF=1
T_rev_1 = zeros((size(sigma,1)),1);

for i = 1:size(sigma,1)
    [T_fwd_0(i,1)]=Forward_Transformation(sigma(i,1),0,P,TP);
    [T_fwd_1(i,1)]=Forward_Transformation(sigma(i,1),1,P,TP);
    [T_rev_0(i,1)]=Reverse_Transformation(sigma(i,1),0,P,TP);
    [T_rev_1(i,1)]=Reverse_Transformation(sigma(i,1),1,P,TP);
end

box on 
plot(T_fwd_0,sigma/(10^6),T_fwd_1,sigma/(10^6),T_rev_0,sigma/(10^6),T_rev_1,sigma/(10^6),'LineWidth',2)
xlabel('Temperature (K)')
ylabel('Stress (MPa)')
% title('SMA Model Phase Diagram')
set(gca,'FontName','Times New Roman','fontsize', 20,'linewidth',1.15)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',3*get(gca,'ticklength'))
legend('\Phi_{fwd, \xi = 0}','\Phi_{fwd, \xi = 1}','\Phi_{rev,  \xi = 0}','\Phi_{rev,  \xi = 1}','Location','southeast')
