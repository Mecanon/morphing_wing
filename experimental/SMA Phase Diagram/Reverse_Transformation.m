function [ T ] = Reverse_Transformation( sigma, MVF,P,TP )
% Function to return the Temperature (T) value for a stress (sigma)under a
% forward transformation at an inputed Martensite Volume Fraction (MVF)

% Calculate the hardening function
f_rev = .5*TP.a2*(1+MVF^P.n3-(1-MVF)^P.n4)-TP.a3;

delta_S=(1/P.E_M-1/P.E_A);


% Output the temperature using the Transformation Surface equation (set
% equal to 0 during transformation)
% Note: assumed MVF_r=1 for Phase Diagram
MVF_r=1;
eps_t_r= H_cursolver( sigma, P.sigma_crit, P.k, P.H_min, P.H_sat);
T = (-(1+TP.D)*sigma*eps_t_r/MVF_r-1/2*delta_S*sigma^2+TP.rho_delta_u0+f_rev-TP.Y_0_t)/TP.rho_delta_s0;


end