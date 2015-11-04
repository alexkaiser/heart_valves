function T = tension_linear(X_current,X_nbr,R_current,R_nbr,k_spr, ref_frac)
% 
% Returns the tension in the linear constitutive law 
% Multiplied by the magnitude of 
% Here we have taken the reciprocal norm in the tangent term
% and multiplied it through 

R_norm = ref_frac * norm(R_current - R_nbr); 
X_norm = norm(X_current - X_nbr); 

T = k_spr * (1/R_norm - 1/X_norm); 
