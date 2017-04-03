function J = build_jacobian_bead_slip(leaflet)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 
    
    X_current                 = leaflet.X; 
    p_0                       = leaflet.p_0; 
    alpha                     = leaflet.alpha; 
    beta                      = leaflet.beta; 
    chordae                   = leaflet.chordae; 
    chordae_idx_left          = leaflet.chordae_idx_left; 
    chordae_idx_right         = leaflet.chordae_idx_right;
    j_max                     = leaflet.j_max; 
    k_max                     = leaflet.k_max; 
    du                        = leaflet.du; 
    dv                        = leaflet.dv; 
    is_internal               = leaflet.is_internal; 
    is_bc                     = leaflet.is_bc; 
    linear_idx_offset         = leaflet.linear_idx_offset; 
    num_trees                 = leaflet.num_trees; 
    total_internal_with_trees = leaflet.total_internal_with_trees; 
    
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

    % there are fewer than 15 nnz per row
    % if using the redundant features on sparse creation use more 
    capacity = 10 * 15 * total_internal_with_trees; 
    
    % build with indices, then add all at once 
    nnz_placed = 0; 
    j_idx      = zeros(capacity, 1); 
    k_idx      = zeros(capacity, 1); 
    vals       = zeros(capacity, 1); 
    
    
    % constant stride arrays for building 3x3 blocks 
    j_offsets = [0 1 2 0 1 2 0 1 2]'; 
    k_offsets = [0 0 0 1 1 1 2 2 2]';

    
    for tree_idx = 1:num_trees
        
        [m N_chordae] = size(chordae(tree_idx).C);         
        free_edge_idx = chordae(tree_idx).free_edge_idx; 
        C             = chordae(tree_idx).C; 
        
        if tree_idx == 1 
            left_side = true; 
        else 
            left_side = false; 
        end 
        
        if left_side
            chordae_idx = chordae_idx_left;  
        else 
            chordae_idx = chordae_idx_right;
        end 
        
        
        % free edge terms first              
        for i=1:size(free_edge_idx, 1)
            j = free_edge_idx(i,1);
            k = free_edge_idx(i,2);

            range_current = linear_idx_offset(j,k) + (1:3); 

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
            
            J_tmp = alpha * dv * tangent_jacobian(X, X_nbr); 
            
            if repulsive_potential
                J_tmp = J_tmp + alpha * dv * du^2 * c_repulsive_circumferential * replusive_jacobian(X,X_nbr,power); 
            end 
            
            if decreasing_tension
                J_tmp = J_tmp + alpha * dv * dec_tension_jacobian(X,X_nbr,du,c_dec_tension_circumferential); 
            end 
            
            % current term is always added in 
            % this gets no sign 
            % this is always at the current,current block in the matrix 
            place_tmp_block(range_current, range_current, J_tmp); 
            
            % If the neighbor is an internal point, it also gets a Jacobian contribution 
            % This takes a sign
            if is_internal(j_nbr,k_nbr)
                range_nbr  = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, -J_tmp); 
            end 

            % interior neighbor is up in k, always  
            j_nbr = j;     
            k_nbr = k+1; 

            % Anterior radial
            X_nbr = X_current(:,j_nbr,k_nbr); 

            J_tmp = beta * du * tangent_jacobian(X, X_nbr); 
            
            if repulsive_potential
                J_tmp = J_tmp + beta * du * dv^2 * c_repulsive_radial * replusive_jacobian(X,X_nbr,power); 
            end 
            
            if decreasing_tension
                J_tmp = J_tmp + beta * du * dec_tension_jacobian(X,X_nbr,dv,c_dec_tension_radial); 
            end

            % current term is always added in 
            % this gets no sign 
            % this is always at the current,current block in the matrix 
            place_tmp_block(range_current, range_current, J_tmp); 

            
            % If the neighbor is an internal point, it also gets a Jacobian contribution 
            % This takes a sign
            if is_internal(j_nbr,k_nbr)
                range_nbr  = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, -J_tmp); 
            end

            % current node has a chordae connection
            if chordae_idx(j,k)

                kappa = chordae(tree_idx).k_0;

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k) + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = C(:,idx_chordae);

                J_tmp = kappa * tangent_jacobian(X, X_nbr); 
                
                if repulsive_potential
                    J_tmp = J_tmp + kappa * c_repulsive_chordae * du^2 * replusive_jacobian(X,X_nbr,power); 
                end 
                
                if decreasing_tension
                    J_tmp = J_tmp + kappa * dec_tension_jacobian(X,X_nbr,du,c_dec_tension_chordae); 
                end

                % current term is always added in 
                % this gets no sign 
                % this is always at the current,current block in the matrix 
                place_tmp_block(range_current, range_current, J_tmp); 
                
                % chordae range 
                range_nbr = chordae(tree_idx).min_global_idx + 3*(idx_chordae-1) + (0:2);
                place_tmp_block(range_current, range_nbr, -J_tmp); 
                
            else
                error('free edge point required to have chordae connection'); 
            end

        end
    
    end 
    



    % Internal anterior leaflet 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max

            % Internal points, not on free edge 
            if is_internal(j,k) && ~chordae_idx_left(j,k) && ~chordae_idx_right(j,k)

                X = X_current(:,j,k); 

                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset(j,k) + (1:3); 

                % pressure portion 
                % always four neighbors
                if p_0 ~= 0.0
                    
                    j_nbr = j+1; 
                    k_nbr = k; 
                    
                    J_pressure = -(p_0/(4*du*dv)) * cross_matrix(X_current(:,j,k+1) - X_current(:,j,k-1)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j-1; 
                    k_nbr = k; 
                    J_pressure =  (p_0/(4*du*dv)) * cross_matrix(X_current(:,j,k+1) - X_current(:,j,k-1)) ; 
                    
                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k+1; 
                    J_pressure =  (p_0/(4*du*dv)) * cross_matrix(X_current(:,j+1,k) - X_current(:,j-1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k-1; 
                    J_pressure = -(p_0/(4*du*dv)) * cross_matrix(X_current(:,j+1,k) - X_current(:,j-1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 

                end 


                for j_nbr = [j-1,j+1]

                    k_nbr = k; 

                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))

                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 

                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 
                        J_tmp = alpha/du * tangent_jacobian(X, X_nbr); 
                        
                        if repulsive_potential
                            J_tmp = J_tmp + alpha/du * c_repulsive_circumferential * du^2 * replusive_jacobian(X,X_nbr,power); 
                        end 
                        
                        if decreasing_tension
                            J_tmp = J_tmp + alpha/du * dec_tension_jacobian(X,X_nbr,du,c_dec_tension_circumferential); 
                        end 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, J_tmp); 

                        % If the neighbor is an internal point, on current leaflet 
                        % it also gets a Jacobian contribution 
                        % This takes a sign
                        if nbr_jacobian_needed 
                            place_tmp_block(range_current, range_nbr, -J_tmp); 
                        end

                    end 
                end

                
                % v tension terms 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 

                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))

                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 

                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                        % There is a 1/dv term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 
                        J_tmp = beta/dv * tangent_jacobian(X, X_nbr); 
                        
                        if repulsive_potential
                            J_tmp = J_tmp + beta/dv * c_repulsive_radial * dv^2 * replusive_jacobian(X,X_nbr,power); 
                        end
                        
                        if decreasing_tension
                            J_tmp = J_tmp + beta/dv * dec_tension_jacobian(X,X_nbr,dv,c_dec_tension_radial); 
                        end
                        
                        % current term is always added in 
                        % this gets no sign  
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, J_tmp); 

                        % If the neighbor is an internal point, on current leaflet 
                        % it also gets a Jacobian contribution 
                        % This takes a sign
                        if nbr_jacobian_needed 
                            place_tmp_block(range_current, range_nbr, -J_tmp); 
                        end 
                        
                    end 
                end
            end
        end
    end

    
    % chordae internal terms         
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 

        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            % this is the same, updating the equations for this component 
            range_current = chordae(tree_idx).min_global_idx + 3*(i-1) + (0:2);               
            % range_chordae(total_internal_leaflet, N_chordae, i, left_side); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

                % if the neighbor is in the chordae 
                if isempty(j_nbr) && isempty(k_nbr) 
                    range_nbr = chordae(tree_idx).min_global_idx + 3*(nbr_idx-1) + (0:2);               
                    % range_chordae(total_internal, N_chordae, nbr_idx, left_side); 
                elseif is_internal(j_nbr, k_nbr)
                    % neighbor is on the leaflet 
                    range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                elseif is_bc(j_nbr, k_nbr)
                    % no block added for neighbor on boundary 
                    range_nbr = []; 
                else 
                    error('Should be impossible, neighbor must be chordae, internal or bc'); 
                end

                % tension Jacobian for this spring 
                J_tmp = k_val * tangent_jacobian(C(:,i), nbr); 
                
                if repulsive_potential
                    J_tmp = J_tmp + k_val * c_repulsive_chordae * du^2 * replusive_jacobian(C(:,i),nbr,power); 
                end
                
                if decreasing_tension
                    J_tmp = J_tmp + k_val * dec_tension_jacobian(C(:,i), nbr, du, c_dec_tension_chordae); 
                end

                % current always gets a contribution from this spring 
                place_tmp_block(range_current, range_current, J_tmp); 

                % range may be empty if papillary muscle, in which case do nothing 
                if ~isempty(range_nbr)
                    place_tmp_block(range_current, range_nbr, -J_tmp); 
                end 
            end 
        end 
    end 
     
 

    J = sparse(j_idx(1:nnz_placed), k_idx(1:nnz_placed), vals(1:nnz_placed), total_internal_with_trees, total_internal_with_trees, nnz_placed);  
    
%%%%%%%%%%
%
% End of real code for building Jacobian 
%  
%%%%%%%%%%


    function place_tmp_block(range_current_loc, range_nbr_loc, block_loc)
        % 
        % Places a 3x3 block in a vectors associated with a sparse matrix 
        % This function is NESTED which means that 
        % it has access to the entire workspace of the calling function 
        % 
        % Input:  
        %   range_current_loc   indices to place in j direction  
        %   range_nbr_loc       indices to place in k direction  
        %   block_loc           block to place 
        % 

        if ~all(size(block_loc) == [3,3])
            error('Must place a 3x3 block'); 
        end 

        % reallocate if too big 
        if nnz_placed + 9 >= capacity
            capacity = 2*capacity; 
            j_idx(capacity) = 0.0; 
            k_idx(capacity) = 0.0; 
            vals(capacity)  = 0.0; 
            fprintf(1, 'Hit reallocation, adjust the initial parameter up so this does not happen.\n'); 
        end 

        j_idx(nnz_placed+1 : nnz_placed+9) = range_current_loc(1) + j_offsets; 
        k_idx(nnz_placed+1 : nnz_placed+9) = range_nbr_loc(1)     + k_offsets; 
        vals (nnz_placed+1 : nnz_placed+9) = block_loc(:); 

        nnz_placed = nnz_placed+9; 

    end 


    function [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor()
        % nested function for getting neighbor 
        % nested functions have full access to current work space 

        if chordae_idx_left(j_nbr,k_nbr) || chordae_idx_right(j_nbr,k_nbr)
            X_nbr = X_current(:,j_nbr,k_nbr); 
            range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
            nbr_jacobian_needed = true; 
        elseif is_internal(j_nbr,k_nbr) 
            X_nbr = X_current(:,j_nbr,k_nbr);
            range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
            nbr_jacobian_needed = true;                        
        elseif is_bc(j_nbr,k_nbr)
            X_nbr = X_current(:,j_nbr,k_nbr); 
            range_nbr = NaN; 
            nbr_jacobian_needed = false; 
        else
            error('requesting nbr in location with no nbr'); 
        end 
    end 

end 

