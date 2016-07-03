function [ MVF, eps_t, E, MVF_r, eps_t_r, eps, Phi_fwd, Phi_rev ] = Implicit_Transformation_Correction_stress( P, TP, chck,MVF, eps_t, E,MVF_r,eps_t_r,sigma, H_cur,eps, T, T_0, Phi_fwd, Phi_rev )
% Function to return the outputs for Martensitic volume ftaction, stress,
% transformation strain, Young's modulus, etc. using implicit integration
% scheme
% Input variables are from elastic prediction

% Return elastic predictions if no transformation is occuring
if chck==0
   
% Forward transformation correction
elseif chck==1
    [ MVF, eps_t, E, MVF_r, eps_t_r, eps, Phi_fwd] = ...
        Implicit_fwd_transformation_correction_stress( MVF, eps_t, E,MVF_r,eps_t_r,sigma, H_cur,eps, T, T_0,...
        Phi_fwd,P,TP );
    
% Reverse transformation correction    
elseif chck==2
    [ MVF, eps_t, E, MVF_r, eps_t_r, eps, Phi_rev] = ...
        Implicit_rev_transformation_correction_stress( MVF, eps_t, E,MVF_r,eps_t_r,sigma, eps, T, T_0,...
        Phi_rev,P, TP);
end
end

