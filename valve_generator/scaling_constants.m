
% Script of hacks on estimating good parameters with bad scaling 
% in form of good scaling. 
%


N = 32; 
n = log2(32); 

du_old = 1.274556735895149e-01; 
dv_old = 1.274556735895149e-01; 

alpha_old = 1; 
beta_old = 1; 

k_0_old = 1.8; 

k_multiplier_old = 1.8; 


k_root_old = 1.889568000000001e+01; 

% Since alpha, beta normalized, 
% alpha, beta, k_0 are 1,1,1.8 


% used directly with this number 
repulsive_coeff_old = 2.238985441466275e-03; 





% New version 

du_new = 1/32;  
dv_new = 1/32; 

alpha_new = 1; 
beta_new  = 1; 

k_m_new = 1.8; 

% leaflet tension coeffs 
% alpha/du 

% connections to free edge on LEAFLET part given 
% alpha * dv
% beta * du 

% same relationship says that k_0 should be mean, multiplied by 1.8 
% and scaled by du and dv 

k_0_32_new = 1.8 * (0.5) * (alpha_new * dv_new + beta_new * du_new)

% leaf base constant
% also equal to the total leaf force in the tree 
k_0_1_new = k_0_32_new * N 

% we know the k multiplier here, can use it to compute root tension 
k_root_new_from_multiplicative_scaling = (k_m_new)^n * k_0_32_new 

% generally k_m is computed to preserve total tension in the root and leaflet 
% apply the formula as a consistency check 
k_m_from_rule = 2.0 * (k_root_new_from_multiplicative_scaling / k_0_1_new)^(1/n)

% repulsive coeff used as 
% alpha/du * coeff * (du^2)
% repulsive_coeff_old = alpha_new/du_new * repulsive_coeff_new * (du_new^2)


repulsive_coeff_new = (du_new/alpha_new) * (1/du_new^2) * repulsive_coeff_old











