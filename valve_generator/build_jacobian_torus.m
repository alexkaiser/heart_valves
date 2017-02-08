function J = build_jacobian_torus(torus)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 
    
    X_current          = torus.X; 
    N                  = torus.N; 
    p_0                = torus.p_0; 
    alpha              = torus.alpha; 
    beta               = torus.beta;  
    du                 = torus.du; 
    dv                 = torus.dv; 
    j_max              = torus.j_max; 
    k_max              = torus.k_max; 
    linear_idx_offset  = torus.linear_idx_offset; 
    
    % total internal points in triangular domain 
    total_internal = 3*N*N; 
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
    

    % Internal leaflet 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max

            X = X_current(:,j,k); 

            % vertical offset does not change while differentiating this equation 
            range_current = linear_idx_offset(j,k) + (1:3); 
            
            % subtract one before taking mod for zero based index 
            j_minus = mod(j-1-1, N) + 1; 
            j_plus  = mod(j+1-1, N) + 1; 
            k_minus = mod(k-1-1, N) + 1; 
            k_plus  = mod(k+1-1, N) + 1;

            % pressure portion 
            % always four neighbors
            if p_0 ~= 0.0

                j_nbr = j_plus; 
                k_nbr = k; 
                J_pressure = -(p_0/(4*du*dv)) * cross_matrix(X_current(:,j,k_plus) - X_current(:,j,k_minus)) ; 

                range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, J_pressure); 
                

                j_nbr = j_minus; 
                k_nbr = k; 
                J_pressure =  (p_0/(4*du*dv)) * cross_matrix(X_current(:,j,k_plus) - X_current(:,j,k_minus)) ; 
                
                range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, J_pressure); 
                

                j_nbr = j; 
                k_nbr = k_plus; 
                J_pressure =  (p_0/(4*du*dv)) * cross_matrix(X_current(:,j_plus,k) - X_current(:,j_minus,k)) ; 

                range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, J_pressure); 
                

                j_nbr = j; 
                k_nbr = k_minus; 
                J_pressure = -(p_0/(4*du*dv)) * cross_matrix(X_current(:,j_plus,k) - X_current(:,j_minus,k)) ; 
                
                range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                place_tmp_block(range_current, range_nbr, J_pressure); 

            end 


            for j_nbr = [j_minus,j_plus]

                k_nbr = k; 

                % X_nbr = X_current(:,j_nbr,k_nbr);
                [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                % There is a 1/du term throughout from taking a finite difference derivative 
                % Place this on the tension variables, one of which apprears in each term 

                tension = alpha/du; 

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

                
            % v tension terms 
            for k_nbr = [k_minus,k_plus]

                j_nbr = j; 

                % X_nbr = X_current(:,j_nbr,k_nbr);
                [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                % There is a 1/du term throughout from taking a finite difference derivative 
                % Place this on the tension variables, one of which apprears in each term 

                tension = beta/dv; 

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

        X_nbr = X_current(:,j_nbr,k_nbr);
        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
        nbr_jacobian_needed = true;                        

    end 

end 

