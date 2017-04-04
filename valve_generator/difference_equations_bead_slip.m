function F = difference_equations_bead_slip(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Requires reference configuration R 
    % 
    % Input
    %     leaflet    Current parameters 
    %
    % Output
    %     F        Values of all difference equation, 3 by triangular array 
    % 

    X_current              = leaflet.X; 
    p_0                    = leaflet.p_0; 
    alpha                  = leaflet.alpha; 
    beta                   = leaflet.beta; 
    chordae                = leaflet.chordae; 
    chordae_idx            = leaflet.chordae_idx; 
    j_max                  = leaflet.j_max; 
    k_max                  = leaflet.k_max; 
    du                     = leaflet.du; 
    is_internal            = leaflet.is_internal; 
    total_internal_leaflet = leaflet.total_internal_leaflet; 
    num_trees              = leaflet.num_trees; 
    
    if num_trees ~= 2
        error('not implemented'); 
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
    

    if isfield(leaflet, 'tension_debug') && leaflet.tension_debug
        tension_debug = true; 
    else 
        tension_debug = false; 
    end 


    F_leaflet = zeros(size(X_current)); 
    
    
    for tree_idx = 1:num_trees

        [m N_chordae] = size(chordae(tree_idx).C);         
        free_edge_idx = chordae(tree_idx).free_edge_idx; 
        C             = chordae(tree_idx).C; 

        if tree_idx == 1 
            left_side = true; 
        else 
            left_side = false; 
        end 

        for i=1:size(free_edge_idx, 1)

            F_tmp = zeros(3,1);

            % left free edge has spring connections up and right on both leaflets
            j = free_edge_idx(i,1);
            k = free_edge_idx(i,2);

            X = X_current(:,j,k); 

            % interior neighbor is right in j on left side, 
            % left in j on right side 
            if left_side
                j_nbr = j + 1;
            else 
                j_nbr = j - 1;
            end 
            k_nbr = k;

            % Anterior circumferential 
            X_nbr = X_current(:,j_nbr,k_nbr); 
            
            % Multiply tension by du to get a force,
            % rather than a force density, here 
            tension = alpha * du; 
            
            if repulsive_potential
                tension = tension - alpha * du * c_repulsive_circumferential * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
            end 
            
            if decreasing_tension
                tension = tension + alpha * du * tension_decreasing(X, X_nbr, du, c_dec_tension_circumferential) ; 
            end 
                
            if tension_debug
                if tension < 0 
                    fprintf('tension = %f, free edge point %d, left = %d, radial\n', tension, i, left_side); 
                end 
            end 
            
            F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            % interior neighbor is up in k, always 
            j_nbr = j;     
            k_nbr = k+1; 

            % Anterior radial
            X_nbr = X_current(:,j_nbr,k_nbr); 
            tension = beta * du; 
            
            if repulsive_potential
                tension = tension - beta * du * c_repulsive_radial * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
            end 
            
            if decreasing_tension
                tension = tension + beta * du * tension_decreasing(X, X_nbr, du, c_dec_tension_radial) ; 
            end
            
            if tension_debug
                if tension < 0 
                    fprintf('tension = %f, free edge point %d, left = %d, circumferential\n', tension, i, left_side); 
                end 
            end 
            
            
            F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            % current node has a chordae connection
            if chordae_idx(j,k).tree_idx == tree_idx 

                kappa = chordae(tree_idx).k_0;

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = C(:,idx_chordae);
                tension = kappa; 
                
                if repulsive_potential
                    tension = tension - kappa * c_repulsive_chordae * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                end 
                
                if decreasing_tension
                    tension = tension + kappa * tension_decreasing(X, X_nbr, du, c_dec_tension_chordae) ; 
                end
                
                if tension_debug
                    if tension < 0 
                        fprintf('tension = %f, free edge point %d, left = %d, chordae\n', tension, i, left_side); 
                    end 
                end 
                
                F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            else
                error('free edge point required to have chordae connection'); 
            end

            F_leaflet(:,j,k) = F_tmp; 

        end 

    end 


    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) && (~chordae_idx(j,k).tree_idx)

                X = X_current(:,j,k); 

                F_tmp = zeros(3,1);

                % pressure term first  
                if p_0 ~= 0
                    F_tmp = F_tmp + (p_0 / 4) * cross(X_current(:,j+1,k) - X_current(:,j-1,k), X_current(:,j,k+1) - X_current(:,j,k-1));                     
                end 

                % u type fibers 
                for j_nbr = [j-1,j+1]

                    k_nbr = k; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 
                    
                    tension = alpha; 
                    
                    if repulsive_potential
                        tension = tension - alpha * c_repulsive_circumferential * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                    end 
                                
                    if decreasing_tension
                        tension = tension + alpha * tension_decreasing(X, X_nbr, du, c_dec_tension_circumferential) ; 
                    end 
                    
                    if tension_debug
                        if tension < 0 
                            fprintf('tension = %f, (j,k) = (%d, %d) radial\n', tension, j, k); 
                        end 
                    end 
                    
                    F_tmp = F_tmp + du * tension * (X_nbr-X)/norm(X_nbr-X); 

                end 

                % v type fibers 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 
                    
                    tension = beta; 
                    
                    if repulsive_potential
                        tension = tension - beta * c_repulsive_radial * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                    end 
                    
                    if decreasing_tension
                        tension = tension + beta * tension_decreasing(X, X_nbr, du, c_dec_tension_radial) ; 
                    end
                    
                    if tension_debug
                        if tension < 0 
                            fprintf('tension = %f, (j,k) = (%d, %d) circumferential\n', tension, j, k); 
                        end 
                    end 
                    
                    
                    F_tmp = F_tmp + du * tension * (X_nbr-X)/norm(X_nbr-X); 

                end 

                F_leaflet(:,j,k) = F_tmp;

            end
        end 
    end
    

    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        F_chordae(tree_idx).vals = zeros(size(C));  

        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 
                
                tension = k_val; 
                
                if repulsive_potential
                    tension = tension - k_val * c_repulsive_chordae * du^2 * power * 1/norm(nbr - C(:,i))^(power+1); 
                end 
                
                if decreasing_tension
                    tension = tension + k_val * tension_decreasing(C(:,i), nbr, du, c_dec_tension_chordae) ; 
                end
                
                if tension_debug
                    if tension < 0 
                        fprintf('tension = %f, (i, nbr_idx, tree_idx) = (%d, %d, %d) chordae\n', tension, i, nbr_idx, tree_idx); 
                    end 
                end 

                tension_by_tangent = tension * (nbr - C(:,i)) / norm(nbr - C(:,i));  

                F_chordae(tree_idx).vals(:,i) = F_chordae(tree_idx).vals(:,i) + tension_by_tangent; 

            end 
        end 
    end 

    F = zeros(total_internal_leaflet, 1); 
    
    % k is required to be the outer loop 
    idx = 1; 
    for k=1:k_max
        for j=1:j_max
            if leaflet.is_internal(j,k)
                F(idx + (0:2)) = F_leaflet(:,j,k); 
                idx = idx + 3; 
            end 
        end 
    end

    for tree_idx = 1:num_trees
        F = [F; F_chordae(tree_idx).vals(:)]; 
    end 
    
end 



