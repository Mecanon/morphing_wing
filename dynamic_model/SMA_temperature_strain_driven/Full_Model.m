function [ sigma, MVF, eps_t, E, MVF_r, eps_t_r, H_cur ] = Full_Model(k, T, eps, P, elastic_check, integration_scheme, MVF_init, eps_t_0, sigma_0, n )
% Function to run the One Dimensional, strain-driven, implicit integration
% scheme

% Current transformation strain at calibration stress
H_cur_cal = H_cursolver(P.sig_cal, P.sig_crit,P.k,P.H_min,P.H_sat);

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
    % H_cur: Current maximum transformational strain
    H_cur = zeros(n,1);
    % eps_t: Transformational Strain
    eps_t = zeros(n,1);
    % sig: Stress
    sigma = zeros(n,1);
    % MVF: Martensitic Volume Fraction
    MVF = zeros(n,1);
    % E: Youngs Modulus
    E = zeros(n,1);
    % eps_t_r: transformational strain at transformation reversal
    eps_t_r = zeros(n,1);
    % MVF_r: Martensic Strain at transformation reversal
    MVF_r = zeros(n,1);
    %Phi_fwd: Forward transformation surface
    Phi_fwd = zeros(n,1);
    %Phi_rev: Reverse transformation surface
    Phi_rev = zeros(n,1);

    % Initialize outputs
    H_cur(1,1) = H_cursolver(sigma_0,P.sig_crit,P.k,P.H_min,P.H_sat);
    eps_t(1,1) = eps_t_0; %P.H_sat; %0;
    sigma(1,1) = sigma_0;
    MVF(1,1) = MVF_init; %0;
    E(1,1)=P.E_M;
else
    load('data.mat',  'sigma','MVF', 'eps_t', 'E', 'MVF_r', 'eps_t_r', 'H_cur', 'Phi_fwd', 'Phi_rev')
    if rem(k,n) == 1
        old_n = size(T,1) - n;
        H_cur = [H_cur(1:old_n,1); zeros(n,1)];
        eps_t = [eps_t(1:old_n,1); zeros(n,1)];
        sigma = [sigma(1:old_n,1); zeros(n,1)];
        MVF = [MVF(1:old_n,1); zeros(n,1)];
        E = [E(1:old_n,1); zeros(n,1)];
        eps_t_r = [eps_t_r(1:old_n,1); zeros(n,1)];
        MVF_r = [MVF_r(1:old_n,1); zeros(n,1)];
        Phi_fwd = [Phi_fwd(1:old_n,1); zeros(n,1)];
        Phi_rev = [Phi_rev(1:old_n,1); zeros(n,1)];
    end
end

% Array for number of iterations required for each load step
% increments = zeros((size(T,1)),1);

%calculate everything for new strain
i = k;

% increments(i,1)=0;
% Initialize Output Variables
MVF(i,1)=MVF(i-1,1);
eps_t(i,1)=eps_t(i-1,1);
E(i,1)=E(i-1,1);
MVF_r(i,1)=MVF_r(i-1,1);
eps_t_r(i,1)=eps_t_r(i-1,1);

% Elastic prediction and transformation Check
% Determine user input for transformation check
if strcmpi(elastic_check,'No') || strcmpi(elastic_check,'N')
    % Non-transformation surface rate-informed selected
    [ sigma(i,1), eps_t(i,1), MVF(i,1), H_cur(i,1), Phi_fwd(i,1), Phi_rev(i,1), chck ] = ...
        Elastic_Transformation_check( P, TP, eps(i,1), T(i,1), T(1,1), E(i,1), MVF(i,1), eps_t(i,1), eps_t_r(i,1), MVF_r(i,1) );

elseif strcmpi(elastic_check,'Yes') || strcmpi(elastic_check,'Y')
    % Transformation surface rate-informed selected
    [ sigma(i,1), eps_t(i,1), MVF(i,1), H_cur(i,1), Phi_fwd(i,1), Phi_rev(i,1), chck ] = ...
        Elastic_Transformation_check_RI(P, TP, eps(i,1), T(i,1), T(1,1), E(i,1), MVF(i,1), eps_t(i,1), eps_t_r(i,1), MVF_r(i,1) );

else
    % Display error if neither Yes or No are selected
    disp('Please input "No" or "Yes" for Elastic Prediction: Transformation Surface Rate-informed');
%     break
end

% Determine User input for transformation correction
if strcmpi(integration_scheme,'explicit') || strcmpi(integration_scheme,'E')
    % Explicit integration scheme selected
    % Display error if n1,n2,n3, or n4 are not equal to 1
    if P.n1 ~= 1 || P.n2 ~=1 || P.n3 ~= 1 || P.n4 ~= 1
        h = msgbox('Smoothness parameters must be changed to 1 for explicit integration scheme.', 'Error','error');
%         break
    end
    % Call function to return output variables for explicit correction
    [ MVF(i,1), eps_t(i,1), E(i,1), MVF_r(i,1), eps_t_r(i,1), sigma(i,1), Phi_fwd(i,1), Phi_rev(i,1) ] = ...
        Explicit_Transformation_Correction(P, TP, chck,...
        MVF(i,1),eps_t(i,1),E(i,1),MVF_r(i,1),eps_t_r(i,1),sigma(i,1), H_cur(i,1),eps(i,1), T(i,1), T(1,1), Phi_fwd(i,1), Phi_rev(i,1) );
elseif strcmpi(integration_scheme,'implicit') || strcmpi(integration_scheme,'I')
    % Implicit integration scheme selected
    % Call function to return output variables for explicit correction
    [ MVF(i,1), eps_t(i,1), E(i,1), MVF_r(i,1), eps_t_r(i,1), sigma(i,1), Phi_fwd(i,1), Phi_rev(i,1) ] = ...
        Implicit_Transformation_Correction(P, TP, chck,...
        MVF(i,1),eps_t(i,1),E(i,1),MVF_r(i,1),eps_t_r(i,1),sigma(i,1), H_cur(i,1),eps(i,1), T(i,1), T(1,1), Phi_fwd(i,1), Phi_rev(i,1) );
else
    % Display error if neither Implicit or Explicit are selected
    disp('Please input "Implicit" or "Explicit" for integration scheme');
end

save('data.mat', 'eps', 'T','sigma','MVF', 'eps_t', 'E', 'MVF_r', 'eps_t_r', 'H_cur', 'Phi_fwd', 'Phi_rev')
end

