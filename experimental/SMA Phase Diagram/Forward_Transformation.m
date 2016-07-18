function [ T ] = Forward_Transformation( sigma, MVF,P,TP )
% Function to return the Temperature (T) value for a stress (sigma)under a
% forward transformation at an inputed Martensite Volume Fraction (MVF)

% Calculate the hardening function
f_fwd = .5*TP.a1*(1+MVF^P.n1-(1-MVF)^P.n2)+TP.a3;

delta_S=(1/P.E_M-1/P.E_A);
H_cur= H_cursolver( sigma, P.sigma_crit, P.k, P.H_min, P.H_sat);

% Output the temperature using the Transformation Surface equation (set
% equal to 0 during transformation)
T = (TP.Y_0_t+f_fwd+TP.rho_delta_u0-1/2*delta_S*sigma^2-(1-TP.D)*abs(sigma)*H_cur)/TP.rho_delta_s0;

end

