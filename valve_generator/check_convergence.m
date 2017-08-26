
path = '/Users/alex/Dropbox/NYU/research/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed'

N_values = 2.^[6:10]; 

order_check(path, N_values)





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






 