function J = build_jacobian_bead_slip(leaflet)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 

    % Copyright (c) 2019, Alexander D. Kaiser
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    X_current                 = leaflet.X; 
    p_0                       = leaflet.p_0; 
    alpha                     = leaflet.alpha; 
    beta                      = leaflet.beta; 
    c_dec_radial              = leaflet.c_dec_radial; 
    c_dec_circumferential     = leaflet.c_dec_circumferential; 
    chordae                   = leaflet.chordae; 
    chordae_idx               = leaflet.chordae_idx; 
    j_max                     = leaflet.j_max; 
    k_max                     = leaflet.k_max; 
    du                        = leaflet.du; 
    is_internal               = leaflet.is_internal; 
    is_bc                     = leaflet.is_bc; 
    linear_idx_offset         = leaflet.linear_idx_offset; 
    num_trees                 = leaflet.num_trees; 
    total_internal_with_trees = leaflet.total_internal_with_trees; 
    
    
    if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
        decreasing_tension = true;  
    else 
        decreasing_tension = false;  
    end

    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
    
    if isfield(leaflet, 'targets_for_bcs') && leaflet.targets_for_bcs 
        targets_for_bcs = true; 
        k_target_net = leaflet.target_net; 
        k_target_papillary = leaflet.target_papillary; 
    else 
        targets_for_bcs = false; 
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


    % Internal anterior leaflet 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max
            
            % Internal points only 
            if is_internal(j,k) 

                X = X_current(:,j,k); 

                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset(j,k) + (1:3); 

                % pressure portion 
                % always four neighbors
                if p_0 ~= 0
                    
                    [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k); 
                    
                    j_nbr = j_plus__1; 
                    k_nbr = k; 
                    
                    J_pressure = -p_0 * m * cross_matrix(X_current(:,j,k_plus__1) - X_current(:,j,k_minus_1)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j_minus_1; 
                    k_nbr = k; 
                    J_pressure =  p_0 * m * cross_matrix(X_current(:,j,k_plus__1) - X_current(:,j,k_minus_1)) ; 
                    
                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k_plus__1; 
                    J_pressure =  p_0 * m * cross_matrix(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 
                    
                    j_nbr = j; 
                    k_nbr = k_minus_1; 
                    J_pressure = -p_0 * m * cross_matrix(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k)) ; 

                    if is_internal(j_nbr,k_nbr)
                        range_nbr = linear_idx_offset(j_nbr,k_nbr) + (1:3);
                        place_tmp_block(range_current, range_nbr, J_pressure); 
                    end 

                end 
                

                for j_nbr_tmp = [j-1,j+1]
                    
                    k_nbr_tmp = k; 
                    
                    [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                    if valid && (~target_k_no_j_spring)

                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 

                        
                        if ~target_spring 
                            alpha_tmp     = alpha(j_spr,k_spr);
                            c_dec_tension = c_dec_circumferential(j_spr,k_spr); 

                            % There is a 1/du term throughout from taking a finite difference derivative 
                            % Place this on the tension variables, one of which apprears in each term 
                            J_tmp = du * alpha_tmp * tangent_jacobian(X, X_nbr); 

                            if decreasing_tension && (alpha_tmp ~= 0)
                                J_tmp = J_tmp + du * alpha_tmp * dec_tension_jacobian(X, X_nbr, du, c_dec_tension); 
                            end 
                        
                        else 
                            error('No j direction targets allowed'); 
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

                
                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 
                    
                    [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid
                        % X_nbr = X_current(:,j_nbr,k_nbr);
                        [X_nbr range_nbr nbr_jacobian_needed] = get_neighbor(); 
                        
                        if ~target_spring 
                            beta_tmp      = beta(j_spr,k_spr); 
                            c_dec_tension = c_dec_radial(j_spr,k_spr); 

                            % There is a 1/du term throughout from taking a finite difference derivative 
                            % Place this on the tension variables, one of which apprears in each term 
                            J_tmp = du * beta_tmp * tangent_jacobian(X, X_nbr); 

                            if decreasing_tension && (beta_tmp ~= 0)
                                J_tmp = J_tmp + du * beta_tmp * dec_tension_jacobian(X, X_nbr, du, c_dec_tension); 
                            end
                            
                        else
                            % connected to a node by a target spring that is a bc 
                            % targets are absolute forces, no du here 
                            J_tmp = tension_zero_rest_length_linear_by_tangent_jacobian(X, X_nbr, k_target_net);                           
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
                
                % current node has a chordae connection
                if chordae_idx(j,k).tree_idx
                    
                    tree_idx = chordae_idx(j,k).tree_idx; 
                    
                    C = chordae(tree_idx).C; 
                    
                    [m N_chordae] = size(chordae(tree_idx).C);
                    c_dec_tension_chordae = chordae(tree_idx).c_dec_chordae_leaf; 
                    du_chordae = 1; 

                    kappa = chordae(tree_idx).k_0;

                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet
                    leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                    % then take the parent index of that number in chordae variables
                    idx_chordae = floor(leaf_idx/2);

                    X_nbr = C(:,idx_chordae);

                    J_tmp = kappa * tangent_jacobian(X, X_nbr); 

                    if decreasing_tension && (kappa ~= 0)
                        J_tmp = J_tmp + kappa * dec_tension_jacobian(X,X_nbr,du_chordae,c_dec_tension_chordae); 
                    end

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
        
        % normalize this, no mesh parameters in chordae computations 
        du_chordae = 1; 
        
        % root is an unknown because it is being treated as a target point 
        % hangle this manually 
        if targets_for_bcs

            root = chordae(tree_idx).root; 
                        
            % root index is zero (and is a separate variable)
            i = 0; 
            
            % root always connects to first point 
            nbr_idx = 1; 
            
            % and current range 
            range_current = range_chordae(chordae, i, tree_idx); 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr R_nbr k_val j_nbr k_nbr c_dec_tension_chordae] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 
        
            % if the neighbor is in the chordae 
            if isempty(j_nbr) && isempty(k_nbr) 
                range_nbr = range_chordae(chordae, nbr_idx, tree_idx); 
            else 
                error('Root neighbor must exist and be in the chordae'); 
            end
            
            % tension Jacobian for this spring 
            J_tmp = k_val * tangent_jacobian(root, nbr); 

            if decreasing_tension && (k_val ~= 0.0)
                J_tmp = J_tmp + k_val * dec_tension_jacobian(root, nbr, du_chordae, c_dec_tension_chordae); 
            end

            % current always gets a contribution from this spring 
            place_tmp_block(range_current, range_current, J_tmp); 

            % range may not be empty here
            place_tmp_block(range_current, range_nbr, -J_tmp); 

            
            % block for attachment to boundary condition with target spring 
            J_tmp = tension_zero_rest_length_linear_by_tangent_jacobian(root, chordae(tree_idx).root_target, k_target_papillary);                             
            place_tmp_block(range_current, range_current, J_tmp); 
            % neighbor block is, by definition, a bc and gets no block
            
        end 
            
        
        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            % this is the same, updating the equations for this component 
            range_current = range_chordae(chordae, i, tree_idx); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val j_nbr k_nbr c_dec_tension_chordae] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

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
                J_tmp = k_val * tangent_jacobian(C(:,i), nbr); 
                
                if decreasing_tension && (k_val ~= 0.0)
                    J_tmp = J_tmp + k_val * dec_tension_jacobian(C(:,i), nbr, du_chordae, c_dec_tension_chordae); 
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
 
        if is_internal(j_nbr,k_nbr) 
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

