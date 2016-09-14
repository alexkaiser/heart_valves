% check_jacobian_taylor_series
% 
% checks jacobian is actually an approximation to the derivative
% to the expected order using a Taylor series 
% 

% reset stream for consistent results 
rand('twister',76599)

epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 

N = 4; 

% Valve skeleton parameters 
valve.r = 1.5; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 


% posterior leaflet data structure 
leaflet.N           = N; 
leaflet.reflect_x   = true; 
leaflet.total_angle = pi; 
leaflet.min_angle   = -leaflet.total_angle/2.0; 
leaflet.max_angle   =  leaflet.total_angle/2.0; 

leaflet.filter.a = 1.0; 
leaflet.filter.h = 2.0; 
leaflet.filter.r = valve.r; 

if leaflet.reflect_x
    leaflet.left_papillary  = [-1; 1; 1] .* valve.left_papillary; 
    leaflet.right_papillary = [-1; 1; 1] .* valve.right_papillary; 
else 
    leaflet.left_papillary  = valve.left_papillary; 
    leaflet.right_papillary = valve.right_papillary; 
end 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = false; 

if ~leaflet.radial_and_circumferential 
    [leaflet.free_edge_idx_left leaflet.free_edge_idx_right] = get_free_edge_ranges(leaflet);
else
    error('Radial and circumferential fibers not implemented ')
end 


% Reference configuration 
[leaflet.R leaflet.is_internal leaflet.is_bc] = build_reference_surface(leaflet); 

% Initial configuration is reference configuration 
leaflet.X = leaflet.R;  

% Spring constants in two directions 
leaflet.alpha    =  1.0; 
leaflet.beta     =  1.0; 
leaflet.p_0      = -2.0; % no pressure for now 
leaflet.ref_frac =  0.5; % generic spring constants reduced by this much 

leaflet.chordae_tree = true; 
if leaflet.chordae_tree
    leaflet.k_0          = 1.0; 
    leaflet.k_multiplier = 2.0; 
    leaflet.tree_frac    = 0.5;
    leaflet.chordae      = add_chordae(leaflet); 
else 
    error('Non-tree chordae not implemented'); 
end 


% eval the difference eqns on the perturbation 
[F F_chordae_left F_chordae_right] = difference_equations(leaflet); 
F_linearized = linearize_internal_points(leaflet, F, F_chordae_left, F_chordae_right);

% jacobian does not change 
J = build_jacobian(leaflet); 

figure; 
spy(J); 
title('jacobian nonzero structure in jacobian tester')

% perturbation also does not change 
Z = zeros(size(leaflet.X)); 
for j=1:N
    for k=1:N
        if leaflet.is_internal(j,k)
            Z(:,j,k) = rand(3,1);  
        end 
    end 
end 

leaflet_Z   = leaflet; 
leaflet_Z.X = Z; 
leaflet_Z.chordae.C_left  = rand(size(leaflet_Z.chordae.C_left)); 
leaflet_Z.chordae.C_right = rand(size(leaflet_Z.chordae.C_right)); 
Z_linearized = linearize_internal_points(leaflet_Z, leaflet_Z.X, leaflet_Z.chordae.C_left, leaflet_Z.chordae.C_right); 



fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    leaflet_perturbation   = leaflet;  
    leaflet_perturbation.X = leaflet.X + ep*leaflet_Z.X; 
    
    leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
    leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

    % eval the difference eqns on the perturbation 
    [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations(leaflet_perturbation); 
    F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

    
    errors(i) = norm(F_perturbed_linearized - F_linearized - ep*J*Z_linearized, 2); 
    
    fprintf('%e\t | %e \n', ep, errors(i)); 

end 

fprintf('\n\n\n\n'); 

figure; 
loglog(epsilon_vals, errors, '-*'); 
hold on 
loglog(epsilon_vals, epsilon_vals.^2, '--'); 

legend('error', 'eps^2')


figure; 

% leaflet part 
for k=1:leaflet.N
    for j=1:leaflet.N
        if leaflet.is_internal(j,k)

            j
            k

            errors = zeros(size(epsilon_vals)); 

            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                leaflet_perturbation   = leaflet;  
                leaflet_perturbation.X = leaflet.X + ep * leaflet_Z.X; 

                leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
                leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

                % eval the difference eqns on the perturbation 
                [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations(leaflet_perturbation); 
                F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = linear_index_offset(j,k,N) + (1:3);                     
                errors(i) = norm(diffs(range)); 

                fprintf('%e\t | %e \n', ep, errors(i)); 

            end 

            fprintf('\n\n\n\n'); 

            loglog(epsilon_vals, errors, '-*'); 
            hold on 
            loglog(epsilon_vals, epsilon_vals.^2, '--'); 

            legend('error', 'eps^2')

        end 
    end 
end 


% chordae part if included 
[m N_chordae] = size(leaflet.chordae.C_left); 
total_internal = 3*N*(N+1)/2; 

for left_side = [true false]
    for i=1:N_chordae

        left_side
        i

        errors = zeros(size(epsilon_vals)); 


        for ep_idx = 1:length(epsilon_vals)

            ep = epsilon_vals(ep_idx); 

            % make a new structure for the perturbation 
            leaflet_perturbation   = leaflet;  
            leaflet_perturbation.X = leaflet.X + ep * leaflet_Z.X; 

            leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
            leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

            % eval the difference eqns on the perturbation 
            [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations(leaflet_perturbation); 
            F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

            diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

            range = range_chordae(total_internal, N_chordae, i, left_side); 
            errors(ep_idx) = norm(diffs(range)); 

            fprintf('%e\t | %e \n', ep, errors(ep_idx)); 

        end 

        fprintf('\n\n\n\n'); 

        loglog(epsilon_vals, errors, '-*'); 
        hold on 
        loglog(epsilon_vals, epsilon_vals.^2, '--'); 

        legend('error', 'eps^2')


    end 
end 
 













