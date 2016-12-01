% check_jacobian_taylor_series
% 
% checks jacobian is actually an approximation to the derivative
% to the expected order using a Taylor series 
% 

% reset stream for consistent results 


N = 8; 

% Initialize structures  
attached = false; 
leaflet_only = true; 
valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only); 

rand('twister',76599)

epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 

leaflet = valve.anterior; 

% eval the difference eqns on the perturbation 
[F_anterior] = difference_equations_bead_slip_leaflet_only(leaflet); 

F_linearized = linearize_internal_points(leaflet, F_anterior);

% jacobian does not change 
J = build_jacobian_bead_slip_leaflet_only(leaflet); 

fig = figure; 
spy(J, 'k'); 
title('Jacobian nonzero structure')

j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 

% perturbation also does not change 
Z = zeros(size(leaflet.X)); 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k)
            Z(:,j,k) = rand(3,1);  
        end 
    end 
end 

leaflet_Z   = leaflet; 
leaflet_Z.X = Z; 
Z_linearized = linearize_internal_points(leaflet_Z, leaflet_Z.X); 


fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    leaflet_perturbation   = leaflet;  
    leaflet_perturbation.X = leaflet.X + ep*leaflet_Z.X; 

    % eval the difference eqns on the perturbation 
    [F_perturbed] = difference_equations_bead_slip_leaflet_only(leaflet_perturbation); 
    F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed); 

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
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)

            j
            k

            errors = zeros(size(epsilon_vals)); 

            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                leaflet_perturbation   = leaflet;  
                leaflet_perturbation.X = leaflet.X + ep * leaflet_Z.X; 

                % eval the difference eqns on the perturbation 
                F_perturbed = difference_equations_bead_slip_leaflet_only(leaflet_perturbation); 
                F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = leaflet_perturbation.linear_idx_offset(j,k) + (1:3);                     
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



