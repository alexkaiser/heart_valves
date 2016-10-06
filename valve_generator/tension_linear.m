function T = tension_linear(X_current,X_nbr,R_current,R_nbr,k_spr,ref_frac)
% 
% Returns the tension in the linear constitutive law 
% 
%

R_norm = ref_frac * norm(R_current - R_nbr); 
X_norm =            norm(X_current - X_nbr); 

T = k_spr * (X_norm/R_norm - 1.0); 
