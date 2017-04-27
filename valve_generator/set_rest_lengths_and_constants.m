function leaflet_with_reference = set_rest_lengths_and_constants(leaflet, strain)
% 
% Assignes spring constants and rest lengths such that the current 
% valve configuration has uniform strain as specified here 
% 
% Tensions are computed from the configuration 
% 


X_current              = leaflet.X; 
alpha                  = leaflet.alpha; 
beta                   = leaflet.beta; 
chordae                = leaflet.chordae; 
chordae_idx            = leaflet.chordae_idx; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
num_trees              = leaflet.num_trees; 

if isfield(leaflet, 'periodic_j')
    periodic_j = leaflet.periodic_j; 
else
    periodic_j = zeros(k_max,1); 
end 


R_u = zeros(j_max, k_max); 
k_u = zeros(j_max, k_max); 
R_v = zeros(j_max, k_max); 
k_v = zeros(j_max, k_max); 


% allocation for tree leaflet connection points 
chordae_with_reference = chordae; 
for tree_idx = 1:num_trees 
    [m N_chordae] = size(chordae_with_reference(tree_idx).C);
    
    % zero spring constants 
    chordae_with_reference(tree_idx).k_vals = zeros(N_chordae,1); 
    
    % allocate rest lengths 
    chordae_with_reference(tree_idx).R_ch   = zeros(N_chordae,1); 
    
    % set parameters that should never be used to empty 
    chordae_with_reference(tree_idx).k_0 = []; 
    
    n_leaves = size(chordae_with_reference(tree_idx).free_edge_idx, 1); 
    
    if n_leaves ~= (N_chordae+1)
        error('Free edge and chordae sizes are inconsistent'); 
    end 
    
    chordae_with_reference(tree_idx).k_free_edge = zeros(n_leaves,1); 
    chordae_with_reference(tree_idx).R_free_edge = zeros(n_leaves,1); 
    
end 


% repulsive potential coefficients, if used 
if isfield(leaflet, 'repulsive_potential') && leaflet.repulsive_potential
    repulsive_potential         = true; 
    power                       = leaflet.repulsive_power; 
    c_repulsive_circumferential = leaflet.c_repulsive_circumferential; 
    c_repulsive_radial          = leaflet.c_repulsive_radial; 
    c_repulsive_chordae         = leaflet.c_repulsive_chordae; 
else 
    repulsive_potential         = false; 
    power                       = 1; 
    c_repulsive_circumferential = 0.0; 
    c_repulsive_radial          = 0.0; 
    c_repulsive_chordae         = 0.0; 
end 

if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
    decreasing_tension = true; 
    c_dec_tension_circumferential = leaflet.c_dec_tension_circumferential; 
    c_dec_tension_radial          = leaflet.c_dec_tension_radial; 
    c_dec_tension_chordae         = leaflet.c_dec_tension_chordae; 
else 
    decreasing_tension = false; 
    c_dec_tension_circumferential = 0.0; 
    c_dec_tension_radial          = 0.0; 
    c_dec_tension_chordae         = 0.0; 
end 
    

% Internal leaflet part 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k)

            X = X_current(:,j,k); 

            % u type fibers 
            % set constants in up direction only here 
            for j_nbr_unreduced = [j-1,j+1] 
                
                % j_spr gets periodic reduction if off the minimum side 
                % meaning it is zero 
                j_spr = min(j, j_nbr_unreduced); 
                if j_spr == 0 
                    j_spr = j_max; 
                end 
                
                % j_nbr may need periodic reduction 
                j_nbr = get_j_nbr(j_nbr_unreduced, k, periodic_j, j_max);
                
                k_nbr = k; 

                k_spr = min(k, k_nbr);
                
                if (j_nbr > 0) && (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
                                
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    tension = alpha; 

                    if repulsive_potential
                        tension = tension - alpha * c_repulsive_circumferential * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                    end 

                    if decreasing_tension
                        tension = tension + alpha * tension_decreasing(X, X_nbr, du, c_dec_tension_circumferential) ; 
                    end 

                    % Here tension is a force per unit length 
                    % So we must multiply by a legnth element to get force 
                    % Take the opposing length element 
                    tension = du * tension; 

                    [k_u(j_spr,k_spr) R_u(j_spr,k_spr)] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 
                
                end 
                
            end 


            % v type fibers 
            for k_nbr = [k-1,k+1] 
                
                % no possible periodicity in j 
                j_nbr = j; 
                
                j_spr = min(j, j_nbr); 
                k_spr = min(k, k_nbr);

                if (j_nbr > 0) && (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
                
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    tension = beta; 

                    if repulsive_potential
                        tension = tension - beta * c_repulsive_radial * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                    end 

                    if decreasing_tension
                        tension = tension + beta * tension_decreasing(X, X_nbr, du, c_dec_tension_radial) ; 
                    end

                    tension = du * tension; 

                    [k_v(j_spr,k_spr) R_v(j_spr,k_spr)] = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 
                
                end 
                
            end 
            
            % current node has a chordae connection
            % data added to free edge array 
            if chordae_idx(j,k).tree_idx
                
                tree_idx = chordae_idx(j,k).tree_idx; 

                [m N_chordae] = size(chordae(tree_idx).C);

                kappa = chordae(tree_idx).k_0;
                
                % index in current free edge array 
                i = chordae_idx(j,k).leaf_idx; 

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = chordae(tree_idx).C(:,idx_chordae);
                tension = kappa; 

                if repulsive_potential
                    tension = tension - kappa * c_repulsive_chordae * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                end 

                if decreasing_tension
                    tension = tension + kappa * tension_decreasing(X, X_nbr, du, c_dec_tension_chordae) ; 
                end
                
                [chordae_with_reference(tree_idx).k_free_edge(i), ... 
                 chordae_with_reference(tree_idx).R_free_edge(i)] ...
                    = get_rest_len_and_spring_constants(X, X_nbr, tension, strain, leaflet); 
                    
            end 
            
        end
    end 
end


% chordae internal terms 
for tree_idx = 1:num_trees

    C = chordae(tree_idx).C; 
    [m N_chordae] = size(C);

    for i=1:N_chordae

        % only go parent wise here
        % child owns parent wise spring constants 
        parent = floor(i/2); 

        nbr_idx = parent;  
        
        % get the neighbors coordinates, reference coordinate and spring constants
        [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

        tension = k_val; 

        if repulsive_potential
            tension = tension - k_val * c_repulsive_chordae * du^2 * power * 1/norm(nbr - C(:,i))^(power+1); 
        end 

        if decreasing_tension
            tension = tension + k_val * tension_decreasing(C(:,i), nbr, du, c_dec_tension_chordae) ; 
        end

        [chordae_with_reference(tree_idx).k_vals(i), ...
         chordae_with_reference(tree_idx).R_ch(i)]  ... 
             = get_rest_len_and_spring_constants(C(:,i), nbr, tension, strain, leaflet); 
   
    end 
end 

% Copy all basic data structures 
leaflet_with_reference = leaflet; 

% Add new information 
leaflet_with_reference.R_u = R_u;
leaflet_with_reference.k_u = k_u;
leaflet_with_reference.R_v = R_v;
leaflet_with_reference.k_v = k_v;

leaflet_with_reference.chordae = chordae_with_reference; 

leaflet_with_reference.diff_eqns = @difference_equations_linear; 
leaflet_with_reference.jacobian  = @build_jacobian_linear; 

leaflet_with_reference.ref_frac = 1.0; 




