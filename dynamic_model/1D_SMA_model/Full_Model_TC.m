function [ sigma, MVF, T, eps_t, E, MVF_r, eps_t_r] = Full_Model_TC(k, t, eps, I, P, elastic_check, T_ambient, sigma_o)
% Function to run the One Dimensional, strain-driven, implicit integration
% scheme
% inputs: time(t), eps(strain), P(material properties), elastic_check(?), T_ambient

% Current transformation strain at calibration stress
H_cur_cal = H_cursolver(P.sig_cal, P.sig_crit,P.k,P.H_sat,P.H_min);

% Partial Derivative of H_cur at calibration stress (dH_cur)
dH_cur=partial_Hcur_sigma(P.sig_cal,P.sig_crit,P.k,P.H_sat,P.H_min);

% Transformation Parameters (structure: TP)
TP.rho_delta_s0 = (-2*(P.C_M*P.C_A)*(H_cur_cal+P.sig_cal*dH_cur+P.sig_cal*(1/P.E_M-1/P.E_A)))/(P.C_M+P.C_A);
TP.D = ((P.C_M-P.C_A)*(H_cur_cal+P.sig_cal*dH_cur+P.sig_cal*(1/P.E_M-1/P.E_A)))/((P.C_M+P.C_A)*(H_cur_cal+...
    P.sig_cal*dH_cur));
TP.a1 = TP.rho_delta_s0*(P.M_f-P.M_s);
TP.a2 = TP.rho_delta_s0*(P.A_s-P.A_f);
TP.a3 = -TP.a1/4*(1+1/(P.n1+1)-1/(P.n2+1))+TP.a2/4*(1+1/(P.n3+1)-1/(P.n4+1));
TP.rho_delta_u0 = TP.rho_delta_s0/2*(P.M_s+P.A_f);
TP.Y_0_t = TP.rho_delta_s0/2*(P.M_s-P.A_f)-TP.a3;

%--------------------------------------------------------------------------
% GENERATION OF MEMORY ARRAYS
%--------------------------------------------------------------------------

% Arrays of output variables
if k == 2
    % T: Temperature
    T = zeros((size(t,1)),1); %cod Pedro não tem temp aqui
    T(1,1)= T_ambient;
    % H_cur: Current maximum transformational strain
    H_cur = zeros((size(t,1)),1);
    % eps_t: Transformational Strain
    eps_t = zeros((size(t,1)),1);
    % sig: Stress
    sigma = zeros((size(t,1)),1);
    % MVF: Martensitic Volume Fraction
    MVF = zeros((size(t,1)),1);
    % E: Youngs Modulus
    E = zeros((size(t,1)),1);
    % eps_t_r: transformational strain at transformation reversal
    eps_t_r = zeros((size(t,1)),1);
    % MVF_r: Martensic Strain at transformation reversa
    MVF_r = zeros((size(t,1)),1);
    %Phi_fwd: Forward transformation surface
    Phi_fwd = zeros((size(t,1)),1);
    %Phi_rev: Reverse transformation surface
    Phi_rev = zeros((size(t,1)),1);
    
    % Initialize outputs
    H_cur(1,1) = H_cursolver(sigma_o, P.sig_crit,P.k,P.H_min,P.H_sat);
    eps_t(1,1) = H_cursolver(sigma_o, P.sig_crit,P.k,P.H_min,P.H_sat);
    sigma(1,1) = sigma_o;
    MVF(1,1) = 1.;
    E(1,1)=P.E_M;
else
     load('data.mat',  'sigma','MVF', 'eps_t', 'E', 'MVF_r', 'eps_t_r', 'H_cur', 'Phi_fwd', 'Phi_rev','T')
end


% Array for number of iterations required for each load step
%increments = zeros((size(t,1)),1);


% Initialize Output Variables
MVF(k,1)= MVF(k-1,1);
eps_t(k,1)=eps_t(k-1,1);
E(k,1)=E(k-1,1);
MVF_r(k,1)=MVF_r(k-1,1);
eps_t_r(k,1)=eps_t_r(k-1,1);
    
% Calculate heat source
rho_E = MVF(k,1)*P.rho_E_M + (1-MVF(k,1))*P.rho_E_A;
r = (rho_E/P.rho)*(4*I(k,1)/(pi*P.d^2))^2;
% Elastic prediction and transformation Check
% Determine user input for transformation check
if strcmpi(elastic_check,'No') || strcmpi(elastic_check,'N')
    % Non-transformation surface rate-informed selected
    [ sigma(k,1), T(k,1), eps_t(k,1), MVF(k,1), H_cur(k,1), Phi_fwd(k,1), Phi_rev(k,1), chck ] = ...
        Elastic_Transformation_check_TC( P, TP, eps(k,1), eps(k-1,1), t(k,1), t(k-1,1), T(k-1,1), T(1,1), sigma(k-1,1), Phi_fwd(k-1,1), Phi_rev(k-1,1), E(k,1), MVF(k,1), eps_t(k,1), eps_t_r(k,1), MVF_r(k,1), r );

elseif strcmpi(elastic_check,'Yes') || strcmpi(elastic_check,'Y')
    % Transformation surface rate-informed selected
    [ sigma(k,1), T(k,1), eps_t(k,1), MVF(k,1), H_cur(k,1), Phi_fwd(k,1), Phi_rev(k,1), chck ] = ...
        Elastic_Transformation_check_RI_TC(P, TP, eps(k,1), eps(k-1,1), t(k,1), t(k-1,1), T(k-1,1), T(1,1), sigma(k-1,1), Phi_fwd(k-1,1), Phi_rev(k-1,1), E(k,1), MVF(k,1), eps_t(k,1), eps_t_r(k,1), MVF_r(k,1), r );

else
    % Display error if neither Yes or No are selected
    disp('Please input "No" or "Yes" for Elastic Prediction: Transformation Surface Rate-informed');
end
    
    % Use inplicit integration scheme considering thermomechanical coupling
    % for transformation correction
    
[ MVF(k,1), T(k,1), eps_t(k,1), E(k,1), MVF_r(k,1), eps_t_r(k,1), sigma(k,1), Phi_fwd(k,1), Phi_rev(k,1) ] = ...
    Implicit_Transformation_Correction_TC(P, TP, chck,...
    MVF(k,1),eps_t(k,1),E(k,1),MVF_r(k,1),eps_t_r(k,1),sigma(k,1), H_cur(k,1),eps(k,1), T(k,1), T(1,1), Phi_fwd(k,1), Phi_rev(k,1) );
   
save('data.mat', 'eps', 'T','sigma','MVF', 'eps_t', 'E', 'MVF_r', 'eps_t_r', 'H_cur', 'Phi_fwd', 'Phi_rev','T')
end

