
N_values = 2.^(6:10); 


% basic, from thesis 
% path = 'meshes/plot_meshes/two_leaflet_8_connector_b7a6aed'

% with enhanced convergence on thesis model 
path = 'meshes/convergence_checks_8_23_18/two_leaflet_8_connector_extra_newton_it_double_prec_917ce58ab'
suffix_name = '_thesis_917ce58ab'; 
N_values = 2.^(5:10); 

% no edge connectors 
% path = 'meshes/plot_meshes/two_leaflet_zero_connector_01e2fb6'
% N_values = 2.^[8:10]; 

% no edge connectors, ratio slightly down 
% path = 'meshes/plot_meshes/two_leaflet_zero_connector_ratio_down_2f5200f'
% suffix_name = '_no_edge_connectors_2f5200f'; 

% path = 'meshes/plot_meshes/two_leaflet_8_connector_new_root_scaling_be020ffc'
% suffix_name = 'new_root_scaling_be020ffc'

% path = 'meshes/plot_meshes/two_leaflet_8_connector_new_root_fixed_ratio_edge_connectors_off_f4f0a94'
% suffix_name = '_new_root_fixed_ratio_edge_connectors_off_f4f0a94'

% path = 'meshes/convergence_checks_8_23_18/two_leaflet_8_connector_new_root_scaling_fixed_ratio_789d69c' 
% suffix_name = '_new_root_scaling_fixed_ratio_789d69c'
% N_values = 2.^[5:11]; 
% 

% path = 'meshes/convergence_checks_8_23_18/cleaf_const_c_root_const_371196f'
% suffix_name = 'cleaf_const_c_root_const_371196f'
% N_values = 2.^[5:10]; 

% path = 'meshes/convergence_checks_8_23_18/cleaf_sqrt_512_over_N_croot_2_on_N_618881d'
% suffix_name = 'cleaf_sqrt_512_over_N_croot_2_on_N_618881d'
% N_values = 2.^(5:11); 

% path = 'meshes/convergence_checks_8_23_18/cleaf_sqrt_N_on_512_16b8f31'; 
% suffix_name = 'cleaf_sqrt_N_on_512_16b8f31'; 
% N_values = 2.^[5:10]; 

% path = 'meshes/convergence_checks_8_23_18/leaf_sqrt_N_on_512_root_2_on_N_cbe9366'
% suffix_name = 'leaf_sqrt_N_on_512_root_2_on_N_cbe9366'
% N_values = 2.^[5:10]; 

% commissural leaflet version, generally terrible 
% path = 'meshes/plot_meshes/comm_leaflets_d599cccee'
% N_values = 2.^[6:9]; 

% initial conditions only to check the check script 
% path = 'initial_meshes'


% path = 'meshes/convergence_checks_8_23_18/cleaf_sqrt_512_over_N_2314694'; 
% suffix_name = 'cleaf_sqrt_512_over_N_2314694'

% path = 'meshes/convergence_checks_8_23_18/cleaf_512_over_N_9d24e2b'
% suffix_name = 'cleaf_512_over_N_9d24e2b'

% path = 'meshes/convergence_checks_8_23_18/cleaf_down_2_b53cd6f'
% suffix_name = 'cleaf_down_2_b53cd6f'

% path = 'meshes/convergence_checks_8_23_18/croot_1_on_n_aaed88b'
% suffix_name = 'croot_1_on_n_aaed88b'

% path = 'meshes/convergence_checks_8_23_18/croot_4_on_n_c6c07c9'
% suffix_name = 'croot_4_on_n_c6c07c9'
% N_values = 2.^(5:11); 


% path = 'meshes/convergence_checks_8_23_18/cleaf_1_on_N_croot_1_on_512_sqrt_N_on_512'; 
% suffix_name = 'cleaf_1_on_N_croot_1_on_512_sqrt_N_on_512'; 
% N_values = 2.^(5:10); 

% path = 'meshes/convergence_checks_8_23_18/cleaf_1_on_N_croot_1_on_512_N_on_512'; 
% suffix_name = 'cleaf_1_on_N_croot_1_on_512_N_on_512'; 
% N_values = 2.^(5:10); 

% path = 'meshes/convergence_checks_8_23_18/const_cleaf_1_on_256_croot_1_on_512'; 
% suffix_name = 'const_cleaf_1_on_256_croot_1_on_512'; 
% N_values = 2.^(5:10); 






order_check(path, N_values, suffix_name)





% radial       = true; 
% bead_slip    = true; 
% attached     = false; 
% leaflet_only = false; 
% optimization = false; 
% repulsive_potential = true; 
% 
% 
% iterations = 6; 
% 
% err_1    = zeros(iterations,1); 
% err_2    = zeros(iterations,1); 
% err_inf  = zeros(iterations,1); 
% 
%  
% N_orig  = 4; 
% N       = N_orig;
% N_vals  = zeros(iterations,1); 
% 
% 
% for it = 1:iterations
%    
%     N_vals(it) = N; 
%     
%     valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential); 
%     
%     p_range = valve.anterior.p_0 .* [0:.1:1]; 
%     
%     [valve pass] = solve_valve(valve, p_range);
%     
%     if ~pass
%         error('Solve failed at N=%d. Exiting convergence test.\n', N); 
%     end 
%     
%     stride = 2; 
%     
%     
%     if it > 1 
%         
%         j_range = 1:stride:valve.anterior.j_max; 
%         
%         % note that k_max is the valve ring and should be discarded 
%         k_range = 2:stride:(valve.anterior.k_max-1); 
%         
%         X_anterior_restriced_j_odd = valve.anterior.X(:, j_range, k_range); 
%         
%         % X_anterior_restriced_j_even = valve.anterior.X(:,1:stride:valve.anterior.j_max, 2:stride:valve.anterior.k_max-1); 
%         
%         diff = abs(X_anterior_restriced_j_odd(:) - X_anterior_prev(:));
%         
%         % replace any nan masks with zeros 
%         diff(isnan(diff)) = 0; 
% 
%         err_1(it)    = norm(diff, 1); 
%         err_2(it)    = norm(diff, 2); 
%         err_inf(it)  = norm(diff, inf);  
%         
%     end 
%     
%     % crop off the boundary condition / valve ring 
%     X_anterior_prev = valve.anterior.X(:,:,1:(valve.anterior.k_max - 1)); 
%     
%     N = N*2; 
% 
% end 
% 
% 
% fprintf('Comparisons of N,2N and 2N,4N points\n'); 
% N = N_orig; 
% for i = 2:it-1
% 
%     order_1    = err_1(i) / err_1(i+1); 
%     order_2    = err_2(i) / err_2(i+1);  
%     order_inf  = err_inf(i) / err_inf(i+1);  
% 
%     fprintf('%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N,  order_1, order_2, order_inf);
% 
% 
%     N = 2*N;
% end
% 






 