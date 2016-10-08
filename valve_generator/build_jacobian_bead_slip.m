function J = build_jacobian_bead_slip(valve)
    % 
    % Builds the Jacobian for the current index and parameter values 
    % 
    % Input 
    %      leaflet   Current parameter values
    % 
    % Output 
    %      J         Jacobian of difference equations 

    anterior  = valve.anterior; 
    posterior = valve.posterior; 
    
    X_anterior        = anterior.X; 
    R_anterior        = anterior.R; 
    p_0_anterior      = anterior.p_0; 
    alpha_anterior    = anterior.alpha; 
    beta_anterior     = anterior.beta; 
    ref_frac_anterior = anterior.ref_frac; 
    C_left            = anterior.chordae.C_left; 
    C_right           = anterior.chordae.C_right; 
    Ref_l             = anterior.chordae.Ref_l; 
    Ref_r             = anterior.chordae.Ref_r; 
    k_0               = anterior.chordae.k_0; 
    chordae_idx_left  = anterior.chordae_idx_left; 
    chordae_idx_right = anterior.chordae_idx_right;
    j_max             = anterior.j_max; 
    k_max             = anterior.k_max; 
    du                = anterior.du; 
    dv                = anterior.dv; 
    is_internal_anterior = anterior.is_internal; 
    is_bc_anterior = anterior.is_bc; 
    free_edge_idx_left   = anterior.free_edge_idx_left; 
    free_edge_idx_right  = anterior.free_edge_idx_right; 
    linear_idx_offset_anterior = anterior.linear_idx_offset; 
    
    X_posterior           = posterior.X; 
    R_posterior           = posterior.R;
    p_0_posterior         = posterior.p_0; 
    alpha_posterior       = posterior.alpha; 
    beta_posterior        = posterior.beta; 
    ref_frac_posterior    = posterior.ref_frac; 
    is_internal_posterior = posterior.is_internal; 
    is_bc_posterior = posterior.is_bc; 
    linear_idx_offset_posterior = posterior.linear_idx_offset; 
    
    [m N_chordae] = size(C_left); 
    
    % total internal points in triangular domain 
    total_internal = 3*(sum(is_internal_anterior(:)) + sum(is_internal_posterior(:))); 
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
    
    % initialize structures for tension variables and jacobians 
    S_anterior_left(k_max-1).val   = 0;  
    S_anterior_left(k_max-1).j     = 0;
    S_anterior_left(k_max-1).k     = 0; 
    S_anterior_left(k_max-1).j_nbr = 0; 
    S_anterior_left(k_max-1).k_nbr = 0; 
    S_anterior_left(k_max-1).G     = zeros(3,1);  
    
    S_anterior_right(k_max-1).val   = 0;  
    S_anterior_right(k_max-1).j     = 0;
    S_anterior_right(k_max-1).k     = 0; 
    S_anterior_right(k_max-1).j_nbr = 0; 
    S_anterior_right(k_max-1).k_nbr = 0; 
    S_anterior_right(k_max-1).G     = zeros(3,1); 
    
    T_anterior(j_max).val   = 0;  
    T_anterior(j_max).j     = 0;
    T_anterior(j_max).k     = 0; 
    T_anterior(j_max).j_nbr = 0; 
    T_anterior(j_max).k_nbr = 0; 
    T_anterior(j_max).G     = zeros(3,1); 
    
    S_posterior_left(k_max-1).val   = 0;  
    S_posterior_left(k_max-1).j     = 0;
    S_posterior_left(k_max-1).k     = 0; 
    S_posterior_left(k_max-1).j_nbr = 0; 
    S_posterior_left(k_max-1).k_nbr = 0; 
    S_posterior_left(k_max-1).G     = zeros(3,1);  
    
    S_posterior_right(k_max-1).val   = 0;  
    S_posterior_right(k_max-1).j     = 0;
    S_posterior_right(k_max-1).k     = 0; 
    S_posterior_right(k_max-1).j_nbr = 0; 
    S_posterior_right(k_max-1).k_nbr = 0; 
    S_posterior_right(k_max-1).G     = zeros(3,1); 
    
    T_posterior(j_max).val   = 0;  
    T_posterior(j_max).j     = 0;
    T_posterior(j_max).k     = 0; 
    T_posterior(j_max).j_nbr = 0; 
    T_posterior(j_max).k_nbr = 0; 
    T_posterior(j_max).G     = zeros(3,1);
                 
    % free edge terms first              
    for i=1:size(free_edge_idx_left, 1)
        j = free_edge_idx_left(i,1);
        k = free_edge_idx_left(i,2);
        
        range_current = linear_idx_offset_anterior(j,k) + (1:3); 
        
        X = X_anterior(:,j,k); 
        R = R_anterior(:,j,k);

        % interior neighbor is right in j 
        j_nbr = j + 1;
        k_nbr = k;

        % Anterior circumferential 
        X_nbr = X_anterior(:,j_nbr,k_nbr); 
        R_nbr = R_anterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_anterior, ref_frac_anterior);

        S_anterior_left(k) = get_linear_tension_struct(X, X_nbr, R, R_nbr, alpha_anterior/du, ref_frac_anterior, j, k, j_nbr, k_nbr); 
        
        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_anterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 
        

        % Posterior circumferential 
        X_nbr = X_posterior(:,j_nbr,k_nbr); 
        R_nbr = R_posterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_posterior, ref_frac_posterior);

        S_posterior_left(k) = get_linear_tension_struct(X, X_nbr, R, R_nbr, alpha_posterior/du, ref_frac_posterior, j, k, j_nbr, k_nbr); 
        
        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_posterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 
        

        % interior neighbor is up in k 
        j_nbr = j;     
        k_nbr = k+1; 

        % Anterior radial
        X_nbr = X_anterior(:,j_nbr,k_nbr); 
        R_nbr = R_anterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, beta_anterior, ref_frac_anterior);
        
        T_anterior(j) = get_linear_tension_struct(X, X_nbr, R, R_nbr, beta_anterior/dv, ref_frac_anterior, j, k, j_nbr, k_nbr); 

        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_anterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 

        % Posterior radial  
        X_nbr = X_posterior(:,j_nbr,k_nbr); 
        R_nbr = R_posterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, beta_posterior, ref_frac_posterior);
        
        T_posterior(j) = get_linear_tension_struct(X, X_nbr, R, R_nbr, beta_posterior/dv, ref_frac_posterior, j, k, j_nbr, k_nbr); 

        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_posterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end

        % current node has a chordae connection
        if chordae_idx_left(j,k)

            kappa = k_0;

            % index that free edge would have if on tree
            % remember that leaves are only in the leaflet
            leaf_idx = chordae_idx_left(j,k) + N_chordae;

            % then take the parent index of that number in chordae variables
            idx_chordae = floor(leaf_idx/2);

            X_nbr = C_left(:,idx_chordae);
            R_nbr = Ref_l (:,idx_chordae);

            J_tension = tension_linear_tangent_jacobian(X,X_nbr,R,R_nbr,kappa,ref_frac_anterior); 

            % current term is always added in 
            % this gets no sign 
            % this is always at the current,current block in the matrix 
            place_tmp_block(range_current, range_current, J_tension); 

            % chordae range 
            left_side = true; 
            range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
            place_tmp_block(range_current, range_nbr, -J_tension); 

        else
            error('free edge point required to have chordae connection'); 
        end
        
    end 
    
    
    % free edge terms first              
    for i=1:size(free_edge_idx_right, 1)
        j = free_edge_idx_right(i,1);
        k = free_edge_idx_right(i,2);
        
        range_current = linear_idx_offset_anterior(j,k) + (1:3); 
        
        X = X_anterior(:,j,k); 
        R = R_anterior(:,j,k);

        % interior neighbor is left in j 
        j_nbr = j - 1;
        k_nbr = k;

        % Anterior circumferential 
        X_nbr = X_anterior(:,j_nbr,k_nbr); 
        R_nbr = R_anterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_anterior, ref_frac_anterior);

        S_anterior_right(k) = get_linear_tension_struct(X, X_nbr, R, R_nbr, alpha_anterior/du, ref_frac_anterior, j, k, j_nbr, k_nbr); 
        
        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_anterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 
        

        % Posterior circumferential 
        X_nbr = X_posterior(:,j_nbr,k_nbr); 
        R_nbr = R_posterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_posterior, ref_frac_posterior);
        
        S_posterior_right(k) = get_linear_tension_struct(X, X_nbr, R, R_nbr, alpha_posterior/du, ref_frac_posterior, j, k, j_nbr, k_nbr); 

        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_posterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 
        

        % interior neighbor is up in k 
        j_nbr = j;     
        k_nbr = k+1; 

        % Anterior radial
        X_nbr = X_anterior(:,j_nbr,k_nbr); 
        R_nbr = R_anterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_anterior, ref_frac_anterior);
        
        T_anterior(j) = get_linear_tension_struct(X, X_nbr, R, R_nbr, beta_anterior/dv, ref_frac_anterior, j, k, j_nbr, k_nbr); 

        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_anterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end 

        % Posterior radial  
        X_nbr = X_posterior(:,j_nbr,k_nbr); 
        R_nbr = R_posterior(:,j_nbr,k_nbr); 
        
        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, alpha_posterior, ref_frac_posterior);

        T_posterior(j) = get_linear_tension_struct(X, X_nbr, R, R_nbr, beta_posterior/dv, ref_frac_posterior, j, k, j_nbr, k_nbr); 
        
        % current term is always added in 
        % this gets no sign 
        % this is always at the current,current block in the matrix 
        place_tmp_block(range_current, range_current, J_tension); 

        % If the neighbor is an internal point, it also gets a Jacobian contribution 
        % This takes a sign
        if is_internal_posterior(j_nbr,k_nbr)
            range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
            place_tmp_block(range_current, range_nbr, -J_tension); 
        end

        % current node has a chordae connection
        if chordae_idx_right(j,k)

            kappa = k_0;

            % index that free edge would have if on tree
            % remember that leaves are only in the leaflet
            leaf_idx = chordae_idx_right(j,k) + N_chordae;

            % then take the parent index of that number in chordae variables
            idx_chordae = floor(leaf_idx/2);

            X_nbr = C_right(:,idx_chordae);
            R_nbr = Ref_r (:,idx_chordae);

            J_tension = tension_linear_tangent_jacobian(X,X_nbr,R,R_nbr,kappa,ref_frac_anterior); 

            % current term is always added in 
            % this gets no sign 
            % this is always at the current,current block in the matrix 
            place_tmp_block(range_current, range_current, J_tension); 

            % chordae range 
            left_side = false; 
            range_nbr = range_chordae(total_internal, N_chordae, idx_chordae, left_side); 
            place_tmp_block(range_current, range_nbr, -J_tension); 

        else
            error('free edge point required to have chordae connection'); 
        end
        
    end 
    
    
    
    % Internal anterior leaflet 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max
            
            % Internal points, not on free edge 
            if is_internal_anterior(j,k) && ~chordae_idx_left(j,k) && ~chordae_idx_right(j,k)

                X = X_anterior(:,j,k); 
                
                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset_anterior(j,k) + (1:3); 


                % pressure portion 
                % zero indexed loop because we are computing indices with mod n 
                if p_0_anterior ~= 0.0
                    
                   error('zero pressure required, non zero not implemented')


%                     % if any index is zero, then 
%                     % the pressure term does not include this triangle
%                     if j_nbr_next && k_nbr_next && j_nbr && k_nbr
% 
%                         % Current has two terms from a product rule 
%                         block     =  - (p_0/6) * cross_matrix(X(:,j_nbr     ,k_nbr     ) - X(:,j,k)) ... 
%                                      + (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));                            
% 
%                         place_tmp_block(range_current, range_current, block); 
% 
% 
%                         % nbr term
%                         % nbr gets differentiated away, and nbr_next stays 
%                         % only added if this is internal 
%                         if is_internal(j_nbr,k_nbr)
%                             range_nbr       = linear_idx_offset(j_nbr,k_nbr) + (1:3);
%                             block = - (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));
%                             place_tmp_block(range_current, range_nbr, block); 
%                         end 
% 
%                         % nbr_next term
%                         % nbr_next gets differentiated away, and nbr stays and gets a sign 
%                         % only added if this is internal 
%                         if is_internal(j_nbr_next,k_nbr_next)
%                             range_nbr_next  = linear_idx_offset(j_nbr_next,k_nbr_next) + (1:3);
%                             block = (p_0/6) * cross_matrix(X(:,j_nbr,k_nbr) - X(:,j,k)); 
%                             place_tmp_block(range_current, range_nbr_next, block); 
%                         end 
%                     end 

                end 


                % u tension terms 
                for j_nbr = [j-1,j+1]
                    
                    k_nbr = k; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal_anterior(j_nbr,k_nbr) || is_bc_anterior(j_nbr,k_nbr))
                    
                        
                        
                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 
                        
                        X_nbr = X_anterior(:,j_nbr,k_nbr);
                        
                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 

                        tension = (1/du) * 0.5 * (S_anterior_left(k).val + S_anterior_right(k).val); 
                        
                        grad_tension_left  = (1/du) * S_anterior_left(k).G; 
                        
                        grad_tension_right = (1/du) * S_anterior_right(k).G; 
                        
                        tangent = (X_nbr - X) / norm(X_nbr - X); 
                        
                        J_tangent = tangent_jacobian(X, X_nbr); 
                        

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        % If the neighbor is an internal point, it also gets a Jacobian contribution 
                        % This takes a sign
                        if is_internal_anterior(j_nbr,k_nbr)
                            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
                            place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
                        end 
                    
                        % Jacobians with respect to inherited tensions
                        J_left = grad_tension_left * tangent'; 

                        j_edge = S_anterior_left(k).j; 
                        k_edge = S_anterior_left(k).k;
                        
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_left); 
                        end 
                        
                        j_nbr_tension = S_anterior_left(k).j_nbr; 
                        k_nbr_tension = S_anterior_left(k).k_nbr;
                        
                        if is_internal_anterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_anterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_left); 
                        end
                        
                        
                        J_right = grad_tension_right * tangent'; 

                        j_edge = S_anterior_right(k).j; 
                        k_edge = S_anterior_right(k).k;
                        
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_right); 
                        end 
                        
                        j_nbr_tension = S_anterior_right(k).j_nbr; 
                        k_nbr_tension = S_anterior_right(k).k_nbr;
                        
                        if is_internal_anterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_anterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_right); 
                        end
                        
                        
                    end 
                end 


                % v tension terms 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal_anterior(j_nbr,k_nbr) || is_bc_anterior(j_nbr,k_nbr))
                        
                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 
                        
                        X_nbr = X_anterior(:,j_nbr,k_nbr);
                        
                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 

                        tension = (1/dv) * T_anterior(j).val; 
                        
                        grad_tension  = (1/dv) * T_anterior(j).G; 
                        
                        tangent = (X_nbr - X) / norm(X_nbr - X); 
                        
                        J_tangent = tangent_jacobian(X, X_nbr); 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        % If the neighbor is an internal point, it also gets a Jacobian contribution 
                        % This takes a sign
                        if (j_nbr > 0) && (k_nbr > 0) && is_internal_anterior(j_nbr,k_nbr)
                            range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
                            place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
                        end
                        
                        
                        % Jacobians with respect to inherited tensions
                        Jac = grad_tension * tangent'; 

                        j_edge = T_anterior(j).j; 
                        k_edge = T_anterior(j).k;
                        
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, Jac); 
                        end 
                        
                        j_nbr_tension = T_anterior(j).j_nbr; 
                        k_nbr_tension = T_anterior(j).k_nbr;
                        
                        if is_internal_anterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_anterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, Jac); 
                        end
                        
                    end 
                end
                
            end
        end
    end
    
    
    % Posterior internal 
    % Note that posterior internal goes to anterior free edge 
    % Zero indices always ignored 
    for j=1:j_max
        for k=1:k_max
            
            % Internal points, not on free edge 
            if is_internal_posterior(j,k) && ~chordae_idx_left(j,k) && ~chordae_idx_right(j,k)

                X = X_posterior(:,j,k); 
                
                % vertical offset does not change while differentiating this equation 
                range_current = linear_idx_offset_posterior(j,k) + (1:3); 


                % pressure portion 
                % zero indexed loop because we are computing indices with mod n 
                if p_0_posterior ~= 0.0
                    
                   error('zero pressure required, non zero not implemented')

                end 


                % u tension terms 
                for j_nbr = [j-1,j+1]
                    
                    k_nbr = k; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal_posterior(j_nbr,k_nbr) || is_bc_posterior(j_nbr,k_nbr))
                        
                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 
  
                        % If point of attachement is a chordae point 
                        % then always on anterior leaflet 
                        % Otherwise, posterior
                        free_edge_nbr = (chordae_idx_left(j_nbr,k_nbr) || chordae_idx_right(j_nbr,k_nbr)); 
                        
                        if free_edge_nbr
                            X_nbr = X_anterior(:,j_nbr,k_nbr);
                        else                             
                            X_nbr = X_posterior(:,j_nbr,k_nbr);
                        end 
                        
                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 

                        tension = (1/du) * 0.5 * (S_posterior_left(k).val + S_posterior_right(k).val); 
                        
                        grad_tension_left  = (1/du) * S_posterior_left(k).G; 
                        
                        grad_tension_right = (1/du) * S_posterior_right(k).G; 
                        
                        tangent = (X_nbr - X) / norm(X_nbr - X); 
                        
                        J_tangent = tangent_jacobian(X, X_nbr); 
                        

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        if free_edge_nbr
                            % connection to free edge on anterior 
                            if is_internal_anterior(j_nbr,k_nbr)
                                range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
                                place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
                            end
                            
                        else                         
                            % connection internal to posterior 
                            if is_internal_posterior(j_nbr,k_nbr)
                                range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
                                place_tmp_block(range_current, range_nbr, -tension * J_tangent); 
                            end 
                        end 
                        
                        
                        % Jacobians with respect to inherited tensions
                        J_left = grad_tension_left * tangent'; 

                        j_edge = S_posterior_left(k).j; 
                        k_edge = S_posterior_left(k).k;
                        
                        % edge term always on anterior 
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_left); 
                        end 
                        
                        % neighbor term always on posterior
                        j_nbr_tension = S_posterior_left(k).j_nbr; 
                        k_nbr_tension = S_posterior_left(k).k_nbr;
                        
                        if is_internal_posterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_posterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_left); 
                        end
                        
                        
                        J_right = grad_tension_right * tangent'; 

                        j_edge = S_posterior_right(k).j; 
                        k_edge = S_posterior_right(k).k;
                        
                        % edge term always on anterior 
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_right); 
                        end 
                        
                        j_nbr_tension = S_posterior_right(k).j_nbr; 
                        k_nbr_tension = S_posterior_right(k).k_nbr;
                        
                        % neighbor term always on posterior
                        if is_internal_posterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_posterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, J_right); 
                        end
                        
                        
                    end 
                end 


                % v tension terms 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    
                    if (j_nbr > 0) && (k_nbr > 0) && (is_internal_posterior(j_nbr,k_nbr) || is_bc_posterior(j_nbr,k_nbr))
                        
                        if chordae_idx_left(j,k) || chordae_idx_right(j,k)
                            error('trying to apply slip model at chordae attachment point'); 
                        end 
                        
                        % If point of attachement is a chordae point 
                        % then always on anterior leaflet 
                        % Otherwise, posterior
                        free_edge_nbr = (chordae_idx_left(j_nbr,k_nbr) || chordae_idx_right(j_nbr,k_nbr)); 
                        
                        if free_edge_nbr
                            X_nbr = X_anterior(:,j_nbr,k_nbr);
                        else                             
                            X_nbr = X_posterior(:,j_nbr,k_nbr);
                        end 
                        
                        
                        
                        % There is a 1/du term throughout from taking a finite difference derivative 
                        % Place this on the tension variables, one of which apprears in each term 

                        tension = (1/dv) * T_posterior(j).val; 
                        
                        grad_tension  = (1/dv) * T_posterior(j).G; 
                        
                        tangent = (X_nbr - X) / norm(X_nbr - X); 
                        
                        J_tangent = tangent_jacobian(X, X_nbr); 

                        % current term is always added in 
                        % this gets no sign 
                        % this is always at the current,current block in the matrix 
                        place_tmp_block(range_current, range_current, tension * J_tangent); 

                        % If the neighbor is an internal point, it also gets a Jacobian contribution 
                        % This takes a sign
                        
                        if free_edge_nbr
                            
                            if (j_nbr > 0) && (k_nbr > 0) && is_internal_anterior(j_nbr,k_nbr)
                                range_nbr  = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
                                place_tmp_block(range_current, range_nbr, -tension * J_tangent);
                            end
                        else
                            if (j_nbr > 0) && (k_nbr > 0) && is_internal_posterior(j_nbr,k_nbr)
                                range_nbr  = linear_idx_offset_posterior(j_nbr,k_nbr) + (1:3);
                                place_tmp_block(range_current, range_nbr, -tension * J_tangent);
                            end
                            
                        end
                        
                        % Jacobians with respect to inherited tensions
                        Jac = grad_tension * tangent'; 

                        j_edge = T_posterior(j).j; 
                        k_edge = T_posterior(j).k;
                        
                        % edge always anterior 
                        if is_internal_anterior(j_edge,k_edge)
                            range_nbr  = linear_idx_offset_anterior(j_edge,k_edge) + (1:3);
                            place_tmp_block(range_current, range_nbr, Jac); 
                        end 
                        
                        j_nbr_tension = T_posterior(j).j_nbr; 
                        k_nbr_tension = T_posterior(j).k_nbr;
                        
                        % internal always posterior
                        if is_internal_posterior(j_nbr_tension,k_nbr_tension)
                            range_nbr  = linear_idx_offset_posterior(k_nbr_tension,k_nbr_tension) + (1:3);
                            place_tmp_block(range_current, range_nbr, Jac); 
                        end
                        
                    end 
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
                [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(anterior, i, nbr_idx, left_side); 

                % if the neighbor is in the chordae 
                if isempty(j_nbr) && isempty(k_nbr) 
                    range_nbr = range_chordae(total_internal, N_chordae, nbr_idx, left_side); 
                elseif is_internal_anterior(j_nbr, k_nbr)
                    % neighbor is on the leaflet 
                    range_nbr = linear_idx_offset_anterior(j_nbr,k_nbr) + (1:3);
                elseif is_bc(j_nbr, k_nbr)
                    % no block added for neighbor on boundary 
                    range_nbr = []; 
                else 
                    error('Should be impossible, neighbor must be chordae, internal or bc'); 
                end

                % tension Jacobian for this spring 
                J_tension = tension_linear_tangent_jacobian(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac_anterior); 

                % current always gets a contribution from this spring 
                place_tmp_block(range_current, range_current, J_tension); 

                % range may be empty if papillary muscle, in which case do nothing 
                if ~isempty(range_nbr)
                    place_tmp_block(range_current, range_nbr, -J_tension); 
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

