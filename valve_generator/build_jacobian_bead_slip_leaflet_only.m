function J = build_jacobian_bead_slip_leaflet_only(leaflet)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 
    
    X_current           = leaflet.X; 
    p_0                 = leaflet.p_0; 
    alpha               = leaflet.alpha; 
    beta                = leaflet.beta; 
    chordae_idx_left    = leaflet.chordae_idx_left; 
    chordae_idx_right   = leaflet.chordae_idx_right;
    j_max               = leaflet.j_max; 
    k_max               = leaflet.k_max; 
    du                  = leaflet.du; 
    dv                  = leaflet.dv; 
    is_internal         = leaflet.is_internal; 
    is_bc               = leaflet.is_bc; 
    linear_idx_offset   = leaflet.linear_idx_offset; 
    
    
    % total internal points in triangular domain 
    total_internal = 3*sum(is_internal(:)); 
    total_points   = total_internal; 

    % there are fewer than 15 nnz per row
    % if using the redundant features on sparse creation use more 
    capacity = 10 * 15 * total_points; 
    
    % build with indices, then add all at once 
    nnz_placed = 0; 
    j_idx      = zeros(capacity, 1); 
    k_idx      = zeros(capacity, 1); 
    vals       = zeros(capacity, 1); 
    
    
    % constant stride arrays for building 3x3 blocks 
    j_offsets = [0 1 2 0 1 2 0 1 2]'; 
    k_offsets = [0 0 0 1 1 1 2 2 2]';
    
    % initialize structures for tension variables and jacobians 
    S_left  = alpha * ones(k_max-1,1); 
    S_right = alpha * ones(k_max-1,1);     
    T       = beta  * ones(j_max,1); 

    % Internal leaflet 
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

                        tension = (1/(2*du)) * (S_left(k) + S_right(k));

                        J_tangent = tangent_jacobian(X, X_nbr); 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        % If the neighbor is an internal point, on current leaflet 
                        % it also gets a Jacobian contribution 
                        % This takes a sign
                        if nbr_jacobian_needed 
                            place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
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

                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 

                        tension = (1/dv) * T(j); 

                        J_tangent = tangent_jacobian(X, X_nbr); 
                        
                        % current term is always added in 
                        % this gets no sign  
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        % If the neighbor is an internal point, on current leaflet 
                        % it also gets a Jacobian contribution 
                        % This takes a sign
                        if nbr_jacobian_needed 
                            place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
                        end 

                    end 
                end
            end
        end
    end
 

    J = sparse(j_idx(1:nnz_placed), k_idx(1:nnz_placed), vals(1:nnz_placed), total_points, total_points, nnz_placed);  
    
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

        if is_bc(j_nbr,k_nbr)
            X_nbr = X_current(:,j_nbr,k_nbr); 
            range_nbr = NaN; 
            nbr_jacobian_needed = false; 
        elseif is_internal(j_nbr,k_nbr) 
            X_nbr = X_current(:,j_nbr,k_nbr);
            range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
            nbr_jacobian_needed = true;                        
        else
            error('requesting nbr in location with no nbr'); 
        end 
    end 

end 

