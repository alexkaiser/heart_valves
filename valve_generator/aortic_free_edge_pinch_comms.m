function leaflet = aortic_free_edge_pinch_comms(leaflet, extra_stretch_radial, N_to_pinch)

    j_max  = leaflet.j_max; 
    k_max  = leaflet.k_max; 
    N_each = leaflet.N_each; 

    R_v    = leaflet.R_v; 
    R_u    = leaflet.R_u; 

    if isfield(leaflet, 'N_leaflets')
        N_leaflets = leaflet.N_leaflets; 
    else 
        N_leaflets = 3; 
    end 
    
    if N_to_pinch < 2 
        error('must pinch at least 2 points for a smooth spline')
    end 
    
    if N_to_pinch > N_each/2
        error('Cannot pinch more than N_each/2 points')
    end 
    
    
    X = leaflet.X; 

    debug = true; 
    debug_text = false; 

    is_bc = leaflet.is_bc;     
    linear_idx_offset         = zeros(j_max, k_max); 
    point_idx_with_bc         = zeros(j_max, k_max); 
    
    free_edge_idx_set = zeros(j_max,1); 
    
    % find free_edge_radius, rest length of half of the free edge 
    free_edge_radius = 0.0; 
    for j_tmp = 1:(N_each/2)
        j_nbr_tmp = j_tmp - 1; 
        k_nbr_tmp = k_max; 
        [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j_tmp, k_max, j_nbr_tmp, k_nbr_tmp); 
        if valid && (~target_spring) && (~target_k_no_j_spring)
            free_edge_radius = free_edge_radius + R_u(j_spr,k_spr); 
        end 
    end

    % find center leaflet height 
    center_leaflet_height = extra_stretch_radial * sum(R_v(N_each/2, 1:(k_max-1))); 

    % annular radius 
    radius = leaflet.r; 

    % height of entire annulus 
    normal_height = leaflet.skeleton.normal_height; 
    
    free_edge_for_initial_conds = sqrt(radius^2 + (normal_height - sqrt(center_leaflet_height^2 - radius^2))^2); 
    
    if free_edge_for_initial_conds > free_edge_radius
        error("requesting positive free edge strain"); 
    end 
    
    for leaflet_idx=1:3

        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (leaflet_idx-1)*N_each;         

        prev_comm_idx = min_idx; 
        if prev_comm_idx == 0
            prev_comm_idx = j_max; 
        end 

        % commissures do not get interpolated 
        free_edge_idx_set(prev_comm_idx) = 1; 
        
        % this point on the leaflet is below the vertical midline and so paired to the previous commissure 
        comm_prev = X(:,prev_comm_idx,k_max); 

        delta_f = free_edge_for_initial_conds / (N_each/2); 

        k=k_max;  

        % attaches free edges 
        for j=1:N_to_pinch

            ring_point = X(:,j + min_idx ,1); 

            % leaflets numbered 1-3
            % comm 1 between leaflet 1,2
            % comm 2 between 2 and 3
            % comm 3 at j_max between leaflet 3 and 1

            % free edge on inside of leaflet 1 only here 

            % leaflet rest height at current point 
            total_height_current = extra_stretch_radial * sum(R_v(j + min_idx, 1:(k_max-1)));

            free_edge_this_j = j * delta_f; 

            % this point on the leaflet is below the vertical midline and so paired to the previous commissure 

            % compute the rest length at the free edge and assure that there's a compressive strain there 
            rest_length_free_edge = 0.0; 
            for j_tmp = (1+min_idx):(j+min_idx)
                j_nbr_tmp = j_tmp - 1; 
                k_nbr_tmp = k_max; 
                [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j_tmp, k_max, j_nbr_tmp, k_nbr_tmp); 
                if valid && (~target_spring) && (~target_k_no_j_spring)
                    rest_length_free_edge = rest_length_free_edge + R_u(j_spr,k_spr); 
                    if (k==k_max) && debug_text
                        fprintf('j_tmp = %d, j_spr = %d, R_u = %f\n', j_tmp, j_spr, R_u(j_spr,k_spr)); 
                    end 
                end 
            end

            if free_edge_this_j > rest_length_free_edge
                error("found current free edge length longer than rest length");             
            end 

            % zero indexed prev_comm_idx minus j, number past the comm 
            j_reflected_temp = mod(prev_comm_idx,j_max) - j; 
            % then set that back with periodicity
            j_reflected = mod(j_reflected_temp,j_max); 
            if j_reflected == 0
                error('this shuold never be zero because zero is the comm point')
            end 

            ring_point_reflected = X(:,j_reflected,1); 

            % leaflet rest height at reflected point
            total_height_reflected = extra_stretch_radial * sum(R_v(j_reflected, 1:(k_max-1))); 

            % relevant distances from each of three points 
            % this is the intersection of three spheres 
            F = @(p) [norm(ring_point - p) - total_height_current; norm(ring_point_reflected - p) - total_height_reflected; norm(comm_prev - p) - free_edge_this_j]; 

            % really drive down that tolerance
            % this is an easy solve 
            options = optimset('Display','off','TolFun',1e-20);
            free_edge_point = fsolve(F,[0;0;0],options);

            if (k==k_max) && debug_text 
                fprintf("residual nonlinear solve for interp point %e = ", norm(F(free_edge_point))); 
            end 

            % based on the rest length 
            X(:,j + min_idx ,k) = free_edge_point; 
            
            % mark this point as set 
            free_edge_idx_set(j + min_idx) = 1; 
            
            % set the corresponding point as well 
            X(:,j_reflected ,k) = free_edge_point; 
            
            % mark this point as set 
            free_edge_idx_set(j_reflected) = 1; 
            
            
            
            if (k==k_max) && debug_text 
                j
                total_height_current
                total_height_reflected
                free_edge_this_j
                free_edge_point
                X_from_interpolation = X(:,j + min_idx ,k)
                fprintf('\n')
            end 
           
        end
        
    end
    
    
    % new free edge length check 
    length_one_leaflet_free_edge = 0; 
    for j=1:N_each
        k=k_max; 
        if ~free_edge_idx_set(j)
            j_nbr_tmp = j+1; 
            k_nbr_tmp = k; 
            [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

            if ~valid 
                error('trying to compute lengths with an invalid rest length')
            end

            length_one_leaflet_free_edge = length_one_leaflet_free_edge + R_u(j_spr,k_spr);                         
        end         
    end 
    fprintf("Circ free edge rest length minus pinched = %f\n", length_one_leaflet_free_edge);
    
    
    
    % set the leaflet center point 
    for comm_idx = 1:N_leaflets

        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (comm_idx-1)*N_each;         

        dj_interp = 1/N_each; 

        prev_comm_idx = min_idx; 
        if prev_comm_idx == 0
            prev_comm_idx = j_max; 
        end 

        comm_prev = X(:,prev_comm_idx,k_max); 
        comm_next = X(:,min_idx + N_each,k_max); 

        j = N_each/2; 
        k = k_max; 
            
        ring_point = X(:,j + min_idx ,1); 

        % total radial rest length of this radial fiber 
        total_rest_length = sum(R_v(j + min_idx, 1:(k-1))); 

        comm_interp_point = (1 - j*dj_interp) * comm_prev ...
                               + j*dj_interp  * comm_next; 

        tangent = (comm_interp_point - ring_point); 
        tangent = tangent / norm(tangent); 

        % based on the rest length 
        X(:,j + min_idx ,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 
        free_edge_idx_set(j + min_idx) = 1; 
    end
    
    
    
    % smooth interpolant to the commissure for other half of leaflet 1 free edge  
    if (N_to_pinch < N_each/2)
        
        for leaflet_idx=1:3

            % point one internal of commissure to point that m
            % N_each is a power of two 
            min_idx = (leaflet_idx-1)*N_each;         

            prev_comm_idx = min_idx; 
            if prev_comm_idx == 0
                prev_comm_idx = j_max; 
            end 

            % min and max internal indices 
            min_j_range = min_idx + 1;
            max_j_range = min_idx + N_each - 1;  
            
            % all the elements for the current leaflet 
            j_range_this_leaflet = [prev_comm_idx, (min_idx+1):(min_idx+N_each)]; 

            % elements that are marked as already set go into the spline 
            % values of j to interpolate over 
            % viewing this as a mapping from index j to value 
            pts = zeros(N_each,1); 
            pts_placed = 0; 
            for j=min_j_range:max_j_range
                if free_edge_idx_set(j)
                    pts_placed = pts_placed + 1; 
                    pts(pts_placed) = j; 
                end 
            end 
            pts = pts(1:pts_placed); 
                        
            for component=1:3
                
                % take the component for 1D interpolation 
                % to reduce headaches 
                vals = X(component,pts,k_max);
                
                for j=min_j_range:max_j_range
                    if ~free_edge_idx_set(j)
                        X(component,j,k) = interp1(pts, vals, j, 'spline'); 
                        % X(component,j,k) = interp1(pts, vals, j); 
                        % free_edge_idx_set(j) = 1; 
                    end 
                end 

            end 
        
        end 
    end 
    
        
    % interpolate internal points of leaflets 
    for leaflet_idx = 1:3
        min_idx = (leaflet_idx-1)*N_each;         

        % dk_interp = 1/k_max; 
        
        for j=1:(N_each-1)
            
            X_ring = X(:,j + min_idx, 1); 
            X_free = X(:,j + min_idx, k_max); 
            
            tangent = (X_free - X_ring); 
            tangent = tangent / norm(tangent); 

            for k=2:k_max
                
                total_rest_length = sum(R_v(j + min_idx, 1:(k-1))); 
                
                % rest length based interpolation 
                X(:,j + min_idx ,k) = total_rest_length * tangent * extra_stretch_radial + X_ring; 
                
                % even interpolation 
                % X(:,j + min_idx, k) = (1 - k*dk_interp) * X_ring + k*dk_interp * X_free; 
            end 
        end         
    end 
    

    dirichlet_free_edge_everywhere = true; 
    if dirichlet_free_edge_everywhere
        for j=1:j_max 
            is_bc(j,k_max) = true; 
        end 
    end 
    
    
    leaflet.X = X; 

    
    if debug 
        figure; 

        x_component = squeeze(X(1,:,:)); 
        y_component = squeeze(X(2,:,:)); 
        z_component = squeeze(X(3,:,:)); 

        width = 1.0; 
        surf(x_component, y_component, z_component, 'LineWidth',width);
        axis equal 
        axis auto 
               
    end 

    % cleanup on utility arrays 
    is_internal = ~is_bc; 

    % Indices for Jacobian building 
    count = 0; 
    for k=1:k_max
        for j=1:j_max
            if is_internal(j,k)
                linear_idx_offset(j,k) = count; 
                count = count + 3; 
            end 
        end 
    end

    % Indices for spring attachments 
    count = 0;
    for k=1:k_max
        for j=1:j_max
            if is_internal(j,k) || is_bc(j,k)
                point_idx_with_bc(j,k) = count; 
                count = count + 1; 
            end 
        end 
    end


    leaflet.is_internal           = is_internal;
    leaflet.is_bc                 = is_bc;
    leaflet.linear_idx_offset     = linear_idx_offset;
    leaflet.point_idx_with_bc     = point_idx_with_bc;

    leaflet.total_internal_leaflet    = 3*sum(leaflet.is_internal(:)); 
    
    
    
end 



