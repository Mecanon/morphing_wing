function [ eps, eps_t, MVF, H_cur, Phi_fwd,Phi_rev, chck ] = Elastic_Transformation_check_stress( P,TP,sigma, eps_prev, T, T_prev, T0, sigma_prev, Phi_fwd_prev, Phi_rev_prev, E, MVF, eps_t, eps_t_r, MVF_r )
% Function to return Elastic prediction and also check for transformation
% without using scaled projections (non-transformation surface rate-informed)

% Returns chck=1 for forward transformation, 2 for reverse transformation 
% and 0 for no transformation

% Solve for strain using elastic prediciton
eps=sigma/E+P.alpha*(T-T0)+eps_t;

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

