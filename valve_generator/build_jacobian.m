function J = build_jacobian(leaflet)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 

    X                 = leaflet.X; 
    R                 = leaflet.R; 
    N                 = leaflet.N; 
    p_0               = leaflet.p_0; 
    alpha             = leaflet.alpha; 
    beta              = leaflet.beta; 
    ref_frac          = leaflet.ref_frac; 
    C_left            = leaflet.chordae.C_left; 
    C_right           = leaflet.chordae.C_right; 
    Ref_l             = leaflet.chordae.Ref_l; 
    Ref_r             = leaflet.chordae.Ref_r; 
    is_internal       = leaflet.is_internal; 
    linear_idx_offset = leaflet.linear_idx_offset; 
    chordae_idx_left  = leaflet.chordae_idx_left; 
    chordae_idx_right = leaflet.chordae_idx_right; 
    
    
    if leaflet.radial_and_circumferential
        error('Radial and circumferential fibers not implemented')
    end 
   
    [m N_chordae] = size(C_left); 
    
    % total internal points in triangular domain 
    total_internal = 3*sum(is_internal(:)); 
    total_points   = total_internal + 3*2*N_chordae; 

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
    
    
    % always 6 pressure neighbors, which may or may not be in bounds
    % relative indices of pressure here 
    % numbered counter clockwise 
    % ignore out of bounds indices, they are not in the pressure 
    % bearing part of the surface 
    pressure_nbrs = [ 0, -1; 
                      1, -1; 
                      1,  0; 
                      0,  1; 
                     -1,  1; 
                     -1,  0]';   

    % All leaflet terms, no chordae 
    % Zero indices always ignored 
    for j=1:N
        for k=1:N
            if leaflet.is_internal(j,k)

                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset(j,k) + (1:3); 


                % pressure portion 
                % zero indexed loop because we are computing indices with mod n 
                if p_0 ~= 0.0
                    for n=0:5

                        j_nbr      = j + pressure_nbrs(1,mod(n  ,6)+1); 
                        k_nbr      = k + pressure_nbrs(2,mod(n  ,6)+1); 
                        j_nbr_next = j + pressure_nbrs(1,mod(n+1,6)+1); 
                        k_nbr_next = k + pressure_nbrs(2,mod(n+1,6)+1);

                        % if any index is zero, then 
                        % the pressure term does not include this triangle
                        if j_nbr_next && k_nbr_next && j_nbr && k_nbr

                            % Current has two terms from a product rule 
                            block     =  - (p_0/6) * cross_matrix(X(:,j_nbr     ,k_nbr     ) - X(:,j,k)) ... 
                                         + (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));                            

                            place_tmp_block(range_current, range_current, block); 


                            % nbr term
                            % nbr gets differentiated away, and nbr_next stays 
                            % only added if this is internal 
                            if is_internal(j_nbr,k_nbr)
                                range_nbr       = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                                block = - (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));
                                place_tmp_block(range_current, range_nbr, block); 
                            end 

                            % nbr_next term
                            % nbr_next gets differentiated away, and nbr stays and gets a sign 
                            % only added if this is internal 
                            if is_internal(j_nbr_next,k_nbr_next)
                                range_nbr_next  = linear_idx_offset(j_nbr_next,k_nbr_next) + (1:3);
                                block = (p_0/6) * cross_matrix(X(:,j_nbr,k_nbr) - X(:,j,k)); 
                                place_tmp_block(range_current, range_nbr_next, block); 
                            end 

                        end
                    end 
                end 


                % u tension terms 
                for j_nbr = [j-1,j+1]
                    
                    k_nbr = k; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (leaflet.is_internal(j_nbr,k_nbr) || leaflet.is_bc(j_nbr,k_nbr))

                        J_tension = tension_jacobian(X(:,j,k),X(:,j_nbr,k_nbr),R(:,j,k),R(:,j_nbr,k_nbr),alpha,ref_frac); 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, J_tension); 

                        % If the neighbor is an internal point, it also gets a Jacobian contribution 
                        % This takes a sign
                        if is_internal(j_nbr,k_nbr)
                            range_nbr  = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                            place_tmp_block(range_current, range_nbr, -J_tension); 
                        end 

% chordae stuff, remove                         
%                     % If neighbor is on the chordae, it has a non zero index 
%                     % This is included here 
%                     elseif idx_chordae ~= 0
%                         range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side);
%                         place_tmp_block(range_current, range_nbr, -J_tension); 

                    
                    end 

                end 


                % v tension terms 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (leaflet.is_internal(j_nbr,k_nbr) || leaflet.is_bc(j_nbr,k_nbr))

                        J_tension = tension_jacobian(X(:,j,k),X(:,j_nbr,k_nbr),R(:,j,k),R(:,j_nbr,k_nbr),beta,ref_frac); 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, J_tension); 

                        % If the neighbor is an internal point, it also gets a Jacobian contribution 
                        % This takes a sign
                        if (j_nbr > 0) && (k_nbr > 0) && is_internal(j_nbr,k_nbr)
                            range_nbr  = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                            place_tmp_block(range_current, range_nbr, -J_tension); 

%                     % If neighbor is on the chordae, it has a non zero index 
%                     % This is included here 
%                     elseif idx_chordae ~= 0
%                         range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
%                         place_tmp_block(range_current, range_nbr, -J_tension); 
                        end 
                    
                    end 

                end
                
                
                if chordae_idx_left(j,k)
                    % current node has a chordae connection 
                    
                    left_side = true; 
                    
                    kappa = leaflet.chordae.k_0; 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_left(j,k) + N_chordae; 
                    
                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  
                    
                    X_nbr = leaflet.chordae.C_left(:,idx_chordae); 
                    R_nbr = leaflet.chordae.Ref_l (:,idx_chordae);
                    
                    J_tension = tension_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,kappa,ref_frac); 
                    
                    % current term is always added in 
                    % this gets no sign 
                    % this is always at the current,current block in the matrix 
                    place_tmp_block(range_current, range_current, J_tension); 
                    
                    % chordae range 
                    range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
                    place_tmp_block(range_current, range_nbr, -J_tension); 
                    
                end 
                
                
                if chordae_idx_right(j,k)
                    % current node has a chordae connection 
                    
                    left_side = false; 
                    
                    kappa = leaflet.chordae.k_0; 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_right(j,k) + N_chordae; 
                    
                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  
                    
                    X_nbr = leaflet.chordae.C_right(:,idx_chordae); 
                    R_nbr = leaflet.chordae.Ref_r  (:,idx_chordae);
                    
                    J_tension = tension_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,kappa,ref_frac); 
                    
                    % current term is always added in 
                    % this gets no sign 
                    % this is always at the current,current block in the matrix 
                    place_tmp_block(range_current, range_current, J_tension); 
                    
                    % chordae range 
                    range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
                    place_tmp_block(range_current, range_nbr, -J_tension); 
                    
                end 

            end
        end
    end

    
    
    
    
    


    % chordae internal terms 
    for left_side = [true false];  

        if left_side
            C = C_left; 
            Ref = Ref_l; 
        else 
            C = C_right; 
            Ref = Ref_r; 
        end 

        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            % this is the same, updating the equations for this component 
            range_current = range_chordae(total_internal, N_chordae, i, left_side); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 

                % if the neighbor is in the chordae 
                if isempty(j_nbr) && isempty(k_nbr) 
                    range_nbr = range_chordae(total_internal, N_chordae, nbr_idx, left_side); 
                else
                    % neighbor is on the leaflet 
                    range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                end 

                % tension Jacobian for this spring 
                J_tension = tension_jacobian(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac); 

                % current always gets a contribution from this spring 
                place_tmp_block(range_current, range_current, J_tension); 

                % range may be empty if papillary muscle, in which case do nothing 
                if ~isempty(range_nbr)
                    place_tmp_block(range_current, range_nbr, - J_tension); 
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

end 

