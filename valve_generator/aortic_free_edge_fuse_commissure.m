function leaflet = aortic_free_edge_fuse_commissure(leaflet, extra_stretch_radial, fused_comm_idx)

    j_max  = leaflet.j_max; 
    k_max  = leaflet.k_max; 
    N_each = leaflet.N_each; 

    R_v    = leaflet.R_v; 
    R_u    = leaflet.R_u; 

    X = leaflet.X; 

    debug = true; 
    debug_text = true; 

    is_bc = leaflet.is_bc; 

    if fused_comm_idx ~= 3
        error('only comm 3 supported for now')
    end 
        
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
    

    leaflet_idx = 1; 

    % point one internal of commissure to point that m
    % N_each is a power of two 
    min_idx = (leaflet_idx-1)*N_each;         

    prev_comm_idx = min_idx; 
    if prev_comm_idx == 0
        prev_comm_idx = j_max; 
    end 

    % this point on the leaflet is below the vertical midline and so paired to the previous commissure 
    comm_prev = X(:,prev_comm_idx,k_max); 

    delta_f = free_edge_for_initial_conds / (N_each/2); 
    
    k=k_max;  
    
    for j=1:(N_each/2)
        
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

        % set to literal zero for symmetry
        tol_component = 1e-14; 
        for component=1:3
            if abs(free_edge_point(component)) < tol_component
                free_edge_point(component) = 0.0; 
            end 
        end 
        
        % based on the rest length 
        X(:,j + min_idx ,k) = free_edge_point; 

        tol_symmetry = 1e-14; 
        if abs(free_edge_point) > tol_symmetry
            error('free edge point not on desired line of symmetry'); 
        end 
        
        
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
    
    % rotate leaflet points for other half of leaflet 1 free edge 
    k=k_max; 
    j_reflected = N_each - 1; 
    theta = 2*pi/3; 
    for j=1:(N_each/2 - 1)
        X(:,j_reflected,k) = rotation_matrix_z(theta) * X(:,j,k);         
        j_reflected = j_reflected - 1;  
    end
    
    % rotate leaflet points for leaflet 3 free edge 
    k=k_max; 
    j_reflected = 2*N_each + 1; 
    theta = 4*pi/3; 
    for j=1:(N_each - 1)
        X(:,j_reflected,k) = rotation_matrix_z(theta) * X(:,j,k);         
        j_reflected = j_reflected + 1;  
    end 
        
    % interpolate remainder of leaflets 
    for leaflet_idx = [1,3]
        min_idx = (leaflet_idx-1)*N_each;         

        dk_interp = 1/k_max; 
        
        for j=2:(N_each-1)
            
            X_ring = X(:,j + min_idx, 1); 
            X_free = X(:,j + min_idx, k_max); 

            for k=2:k_max
                X(:,j + min_idx, k) = (1 - k*dk_interp) * X_ring + k*dk_interp * X_free; 
            end 
        end         
    end 
    
    
    
    % error checking        
    tol_err_check = 1e-14; 
    k=k_max; 
    for leaflet_idx = 1:3            

        min_idx = (leaflet_idx-1)*N_each;     

        prev_comm_idx = min_idx; 
        if prev_comm_idx == 0
            prev_comm_idx = j_max; 
        end 

        for j=1:(N_each/2)
            if (((fused_comm_idx == 1) && (leaflet_idx == 2)) || ...
                ((fused_comm_idx == 2) && (leaflet_idx == 3)) || ...
                ((fused_comm_idx == 3) && (leaflet_idx == 1)))

                % setting 
                j_reflected_temp = mod(prev_comm_idx,j_max) - j; 
                j_reflected = mod(j_reflected_temp,j_max); 

                if norm(X(:,j + min_idx,k_max) - X(:,j_reflected,k_max)) > tol_err_check
                    pt_idx_j = j + min_idx
                    point = X(:,j + min_idx,k_max)
                    nbr_idx_j = j_reflected 
                    neighbor = X(:,j_reflected,k_max)
                    difference = norm(X(:,j + min_idx,k_max) - X(:,j_reflected,k_max))
                    warning('fused points differ by more than tolerance')
                else 
                    % assign to be literally floating point equal 
                    X(:,j_reflected,k_max) = X(:,j + min_idx,k_max); 
                end 

            end                 
        end 
    end                 
     
    

    leaflet.X = X; 

    for j=1:j_max 
        if ~is_bc(j,k)
            error('did not set all bcs')
        end 
    end 

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

end 



