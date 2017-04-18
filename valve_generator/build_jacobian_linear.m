function J = build_jacobian_linear(leaflet)
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
    chordae                   = leaflet.chordae; 
    chordae_idx               = leaflet.chordae_idx; 
    j_max                     = leaflet.j_max; 
    k_max                     = leaflet.k_max;  
    is_internal               = leaflet.is_internal; 
    is_bc                     = leaflet.is_bc; 
    linear_idx_offset         = leaflet.linear_idx_offset; 
    num_trees                 = leaflet.num_trees; 
    total_internal_with_trees = leaflet.total_internal_with_trees;  
    R_u                       = leaflet.R_u;
    k_u                       = leaflet.k_u;
    R_v                       = leaflet.R_v;
    k_v                       = leaflet.k_v;

    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
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

    


%             % current node has a chordae connection
%             if chordae_idx(j,k)
% 
%                 % index that free edge would have if on tree
%                 % remember that leaves are only in the leaflet
%                 leaf_idx = chordae_idx(j,k) + N_chordae;
% 
%                 % then take the parent index of that number in chordae variables
%                 idx_chordae = floor(leaf_idx/2);
% 
%                 X_nbr = C(:,idx_chordae);
% 
%                 J_tmp = tension_linear_tangent_jacobian(X,X_nbr,R_free_edge(i),k_free_edge(i));
% 
%                 % current term is always added in 
%                 % this gets no sign 
%                 % this is always at the current,current block in the matrix 
%                 place_tmp_block(range_current, range_current, J_tmp); 
%                 
%                 % chordae range 
%                 range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
%                 place_tmp_block(range_current, range_nbr, -J_tmp); 
%                 
%             else
%                 error('free edge point required to have chordae connection'); 
%             end



    % Internal anterior leaflet 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max

            % Internal points
            if is_internal(j,k)

                X = X_current(:,j,k); 

                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset(j,k) + (1:3); 

                % pressure portion 
                % always four neighbors
                if (~is_bc(j,k)) && (~chordae_idx(j,k).tree_idx) && (p_0 ~= 0)
                    
                    % periodic reduction of nbr indices 
                    j_plus__1 = get_j_nbr(j+1, k, periodic_j, j_max); 
                    j_minus_1 = get_j_nbr(j-1, k, periodic_j, j_max);
                    
                    
                    j_nbr = j_plus__1; 
                    k_nbr = k; 
                    
                    J_pressure = -(p_0/4) * cross_matrix(X_current(:,j,k+1) - X_current(:,j,k-1)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j_minus_1; 
                    k_nbr = k; 
                    J_pressure =  (p_0/4) * cross_matrix(X_current(:,j,k+1) - X_current(:,j,k-1)) ; 
                    
                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k+1; 
                    J_pressure =  (p_0/4) * cross_matrix(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k-1; 
                    J_pressure = -(p_0/4) * cross_matrix(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 

                end 


                for j_nbr_unreduced = [j-1,j+1]
                    
                    % j_spr gets periodic reduction if off the minimum side 
                    % meaning it is zero 
                    j_spr = min(j, j_nbr_unreduced); 
                    if j_spr == 0 
                        j_spr = j_max; 
                    end 
                    
                    j_nbr = get_j_nbr(j_nbr_unreduced, k, periodic_j, j_max); 

                    k_nbr = k; 
                    k_spr = min(k, k_nbr);

                    if (j_nbr > 0) && (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
                        
                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                        J_tmp = tension_linear_tangent_jacobian(X,X_nbr,R_u(j_spr,k_spr),k_u(j_spr,k_spr));                    

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
                    
                    j_spr = min(j, j_nbr); 
                    k_spr = min(k, k_nbr);

                    if (j_nbr > 0) && (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
                        
                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 
 
                        J_tmp = tension_linear_tangent_jacobian(X,X_nbr,R_v(j_spr,k_spr),k_v(j_spr,k_spr));
                        
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
                
                % current node has a chordae connection
                if chordae_idx(j,k).tree_idx
                    
                    tree_idx = chordae_idx(j,k).tree_idx; 
                    
                    C = chordae(tree_idx).C; 
                    
                    [m N_chordae] = size(chordae(tree_idx).C);
                    
                    % index in current free edge array 
                    i = chordae_idx(j,k).leaf_idx;

                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet
                    leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                    % then take the parent index of that number in chordae variables
                    idx_chordae = floor(leaf_idx/2);

                    X_nbr = C(:,idx_chordae);

                    J_tmp = tension_linear_tangent_jacobian(X,X_nbr,chordae(tree_idx).R_free_edge(i),chordae(tree_idx).k_free_edge(i));

                    % current term is always added in 
                    % this gets no sign 
                    % this is always at the current,current block in the matrix 
                    place_tmp_block(range_current, range_current, J_tmp); 

                    % chordae range 
                    range_nbr = range_chordae(chordae, idx_chordae, tree_idx); 
                    place_tmp_block(range_current, range_nbr, -J_tmp); 
                    
                end 
                    
            end
        end
    end

    
    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        [m N_chordae] = size(C);

        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            % this is the same, updating the equations for this component 
            range_current = range_chordae(chordae, i, tree_idx); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

                % if the neighbor is in the chordae 
                if isempty(j_nbr) && isempty(k_nbr) 
                    range_nbr = range_chordae(chordae, nbr_idx, tree_idx); 
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
                J_tmp = tension_linear_tangent_jacobian(C(:,i),nbr,R_nbr,k_val); 

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

        if chordae_idx(j_nbr,k_nbr).tree_idx
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

