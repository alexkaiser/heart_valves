function [] = rewrite_files_from_solved_valve(N)

base_name = 'mitral_tree'; 
base_name = strcat(base_name, sprintf('_%d', N)); 

load(strcat(base_name, '_final_data')); 

L = 2.5;

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
ratio = 1.5; 


% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
refinement = N/32.0; 


p_physical = 100; 

target_multiplier = 40; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
n_lagrangian_tracers = 8; 


X_config_is_reference = true; 

output_to_ibamr_format(base_name, L, ratio, params_posterior, filter_params_posterior, params_anterior, p_physical, target_multiplier, refinement, n_lagrangian_tracers, X_config_is_reference); 
