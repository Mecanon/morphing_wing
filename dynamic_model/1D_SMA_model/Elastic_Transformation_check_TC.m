function [ sigma, T, eps_t, MVF, H_cur, Phi_fwd,Phi_rev, chck ] = Elastic_Transformation_check_TC( P,TP,eps, eps_prev, t, t_prev, T_prev, T0, sigma_prev, Phi_fwd_prev, Phi_rev_prev, E, MVF, eps_t, eps_t_r, MVF_r, r )
% Function to return Elastic prediction considering thermomechanical 
% couplingand also check for transformation without using scaled 
% projections (non-transformation surface rate-informed)

% Returns chck=1 for forward transformation, 2 for reverse transformation 
% and 0 for no transformation

% Predict change in stress and temperature using system of equations
% obtained from linearization of stress and conservation of energy
% equations
if P.h == 1
    % Gravity:
    g = 9.8; %ms-2
    % Atmospheric pressure
    P_air = 101325.; % Pa
    % Molar
    M = 0.0289644 ; %kg/mol
    % Ideal gas constant
    R = 8.31447 ; %J/(mol K)
    % Air density:
    rho_air = P_air*M / (R*T0);
    % Sutherland's law coefficients
    C1 = 1.458e-6; %kg/m.s.sqrt(K)
    C2 = 110.4; %K
    % Air dynamic viscosity:
    mu_air = (C1 * T0^(3./2)) / (T0+C2);
    % Air kinematic viscosity:
    nu_air = mu_air/rho_air;
    % Air specific heat at constant pressure
    Cp_air = 1.005;
    % Air conductivity
    k_air = 0.0264;
    % Nusselt number coefficients
    alpha_1 = 1.;
    alpha_2 = 0.287;
    % Grashof number for external flow around a cylinder
    Gr = 2*abs(T_prev - T0)/(T_prev + T0)*(g*P.d^3)/(nu_air^2);
    % Prandtl number definition
    Pr = mu_air*Cp_air/k_air;
    % Nusselt number and parameter
    Nu = (alpha_1 + alpha_2*(Gr*Pr/(1 + (0.56/Pr)^(9./16))^(16./9))^(1./6))^2;
    % Calculate convection coefficient h from definition of Nusselt number
    h = k_air*Nu/P.d;
else
    h = 0;
end
% q: Body heat source/sink (for circular cross section)
q=-4*h/P.d*(T_prev-T0);
% Evaluate the change in strain and time over the step
delta_eps=eps-eps_prev;
delta_t=t-t_prev;
% Define matrix 'A' representing 2x2 matrix in linearized stress and energy
% equations
A= [1, E*P.alpha;
    T_prev*P.alpha, P.rho*P.c];
% Evaluate change in stress and temperature ('C' is a dummy vector for the
% change in stress and temperature)
C= inv(A)*[E*delta_eps; (q+P.rho*r)*delta_t];
delta_sigma=C(1,1);
delta_T=C(2,1);
% Use the increments in stress and temperature to calculate the current
% value of each
sigma=sigma_prev+delta_sigma;
T=T_prev+delta_T;

%  Solve for current transformational strain
H_cur=H_cursolver(sigma,P.sig_crit,P.k,P.H_min,P.H_sat);

% Solve for hardening functions 
f_fwd = .5*TP.a1*(1+MVF^P.n1-(1-MVF)^P.n2)+TP.a3;
f_rev = .5*TP.a2*(1+MVF^P.n3-(1-MVF)^P.n4)-TP.a3;

% Solve for the forward transformation surface
chck=0;
Phi_fwd = (1-TP.D)*abs(sigma)*H_cur+.5*(1/P.E_M-1/P.E_A)*sigma^2+...
    TP.rho_delta_s0*T-TP.rho_delta_u0-f_fwd-TP.Y_0_t;
% If the forward transformation surface is greater than zero the material
% is undergoing a forward transformation (chck=1)
if Phi_fwd >0
    chck=1;
end

% Solve for the reverse transformation surface
Phi_rev=0;
if MVF_r == 0
    % Reverse transformation strain remains zero if Martensitic Strain at
    % transformation reversal is zero
else
    Phi_rev = -(1+TP.D)*sigma*eps_t_r/MVF_r-.5*(1/P.E_M-1/P.E_A)...
        *sigma^2-TP.rho_delta_s0*T+TP.rho_delta_u0+f_rev-TP.Y_0_t;
    % If the revverse transformation surface is greater than zero the material
    % is undergoing a forward transformation (chck=2)
    if Phi_rev > 0
        chck=2;
    end
end


end

