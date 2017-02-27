function T = tension_linear(X_current,X_nbr,R,k_spr,ref_frac)
% 
% Returns the tension in the linear constitutive law 
% 
%

if ~exist('ref_frac', 'var')
    ref_frac = 1; 
end 

R      = ref_frac * R; 
X_norm = norm(X_current - X_nbr); 

T = k_spr * (X_norm/R - 1.0); 
