

file_name = 'meshes/aortic_semifinal_6fa5a93/aortic_384_final_data.mat'; 
load(file_name, 'valve', 'valve_with_reference', 'N')


E_circ_star = .15; 
E_rad_star  = .54; 

% exponential rates 
b_circ = 57.456509400487398; 
b_rad = 22.397200094241359; 

tension_circ_no_kappa = exp(b_circ * E_circ_star) - 1.0 

tension_rad_no_kappa = exp(b_rad * E_rad_star) - 1.0 

kappa_circ = valve_with_reference.leaflets.k_u; 

kappa_rad = valve_with_reference.leaflets.k_v; 

non_zero_kappa_circ_idx = (kappa_circ ~= 0); 

non_zero_kappa_rad_idx = (kappa_rad ~= 0); 

kappa_circ_nonzero = kappa_circ(non_zero_kappa_circ_idx); 

min_kappa_circ = min(min(kappa_circ_nonzero))
max_kappa_circ = max(max(kappa_circ_nonzero))
mean_kappa_circ = mean(kappa_circ_nonzero(:))

min_circ_tension = min_kappa_circ * tension_circ_no_kappa
max_circ_tension = max_kappa_circ * tension_circ_no_kappa
mean_circ_tension = mean_kappa_circ * tension_circ_no_kappa

kappa_rad_nonzero = kappa_rad(non_zero_kappa_rad_idx); 


% [min_kappa_circ_index_i, min_kappa_circ_index_j]  = find(min_kappa_circ == kappa_circ)
% [max_kappa_circ_index_i, max_kappa_circ_index_j]  = find(max_kappa_circ == kappa_circ)

min_kappa_rad = min(min(kappa_rad_nonzero))
max_kappa_rad = max(max(kappa_rad_nonzero))
mean_kappa_rad = mean(kappa_rad_nonzero(:))


min_rad_tension = min_kappa_rad * tension_rad_no_kappa
max_rad_tension = max_kappa_rad * tension_rad_no_kappa
mean_rad_tension = mean_kappa_rad * tension_rad_no_kappa

% [min_kappa_rad_index_i, min_kappa_rad_index_j] = find(min_kappa_rad == kappa_rad)
% [max_kappa_rad_index_i, max_kappa_rad_index_j] = find(max_kappa_rad == kappa_rad)
