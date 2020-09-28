function leaflet = aortic_free_edge_fuse_commissure(leaflet, extra_stretch_radial, fused_comm_idx)

    j_max  = leaflet.j_max; 
    k_max  = leaflet.k_max; 
    N_each = leaflet.N_each; 

    R_v    = leaflet.R_v; 
    R_u    = leaflet.R_u; 

    full_leaflet_interp = true; 
    if full_leaflet_interp
        k_range = 2:k_max; 
    else
        k_range = k_max; 
    end 

    X = leaflet.X; 

    debug = true; 
    debug_text = true; 

    is_bc = leaflet.is_bc; 


        
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
    center_leaflet_height = sum(R_v(N_each/2, 1:(k_max-1))) 

    % annular radius 
    radius = leaflet.r

    % height of entire annulus 
    normal_height = leaflet.skeleton.normal_height

    radial_stretch_center = (1/center_leaflet_height) * sqrt(radius^2 + (normal_height - sqrt(free_edge_radius^2 - radius^2))^2) 

    if radial_stretch_center < 1
        error('found compressive radial stretch')
    end 
    if radial_stretch_center > extra_stretch_radial
        error('required radial stretch larger than provided stretch for other leaflets')
    end 

    center_point_from_formula = [0,0,normal_height - sqrt(free_edge_radius^2 - radius^2)]


    for leaflet_idx = 1:3

        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (leaflet_idx-1)*N_each;         

        dj_interp = 1/N_each; 

        prev_comm_idx = min_idx; 
        if prev_comm_idx == 0
            prev_comm_idx = j_max; 
        end 
        
        next_comm_idx = min_idx + N_each; 

        comm_prev = X(:,prev_comm_idx,k_max); 
        comm_next = X(:,min_idx + N_each,k_max); 

        for j=1:(N_each-1)
            for k=k_range
                
                ring_point = X(:,j + min_idx ,1); 
                
                if (((fused_comm_idx == 1) && (leaflet_idx == 1)) || ...
                    ((fused_comm_idx == 1) && (leaflet_idx == 2)) || ...
                    ((fused_comm_idx == 2) && (leaflet_idx == 2)) || ...
                    ((fused_comm_idx == 2) && (leaflet_idx == 3)) || ...
                    ((fused_comm_idx == 3) && (leaflet_idx == 1)) || ...
                    ((fused_comm_idx == 3) && (leaflet_idx == 3)))

                    % leaflets numbered 1-3
                    % comm 1 between leaflet 1,2
                    % comm 2 between 2 and 3
                    % comm 3 at j_max between leaflet 3 and 1
                    
                    % point is total_rest_length * extra_stretch_radial from the ring
                    % same distance from the mirrored point across the closest commissure                     
                    % and rest length along the free edge from the commissure 
                    
                    % leaflet rest height at current point 
                    radial_stretch_j = radial_stretch_local(j, N_each, radial_stretch_center, extra_stretch_radial);
                    total_height_current = radial_stretch_j * sum(R_v(j + min_idx, 1:(k_max-1))); 
                    
                    % this point on the leaflet is below the vertical midline and so paired to the previous commissure 
                    if j <= (N_each/2)
                        comm_for_fused_edge = comm_prev; 
                        
                        total_rest_length_free_edge = 0.0; 
                        
                        for j_tmp = (1+min_idx):(j+min_idx)
                            j_nbr_tmp = j_tmp - 1; 
                            k_nbr_tmp = k_max; 
                            [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j_tmp, k_max, j_nbr_tmp, k_nbr_tmp); 
                            if valid && (~target_spring) && (~target_k_no_j_spring)
                                total_rest_length_free_edge = total_rest_length_free_edge + R_u(j_spr,k_spr); 
                                if (k==k_max) && debug_text
                                    fprintf('j_tmp = %d, j_spr = %d, R_u = %f\n', j_tmp, j_spr, R_u(j_spr,k_spr)); 
                                end 
                            end 
                        end
                        
                        % total_rest_length_free_edge = sum(R_v( min_idx:(j+min_idx), k_max)); 
                
                        % zero indexed prev_comm_idx minus j, number past the comm 
                        j_reflected_temp = mod(prev_comm_idx,j_max) - j; 
                        
                        % then set that back with periodicity
                        j_reflected = mod(j_reflected_temp,j_max); 
                        
                        if j_reflected == 0
                            error('this shuold never be zero because zero is the comm point')
                        end 
                        
                        ring_point_reflected = X(:,j_reflected,1); 
                        
                        % leaflet rest height at reflected point
                        radial_stretch_j_reflected = radial_stretch_local(j_reflected, N_each, radial_stretch_center, extra_stretch_radial);
                        total_height_reflected = radial_stretch_j_reflected * sum(R_v(j_reflected, 1:(k_max-1))); 
                                                                        
                    else 
                        comm_for_fused_edge = comm_next; 
                        
                        total_rest_length_free_edge = 0.0; 
                        
                        for j_tmp = (j+min_idx):(N_each+min_idx-1)
                            j_nbr_tmp = j_tmp + 1; 
                            k_nbr_tmp = k_max; 
                            [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j_tmp, k_max, j_nbr_tmp, k_nbr_tmp); 
                            if valid && (~target_spring) && (~target_k_no_j_spring)
                                total_rest_length_free_edge = total_rest_length_free_edge + R_u(j_spr,k_spr); 
                                if (k==k_max) && debug_text
                                    fprintf('j_tmp = %d, j_spr = %d, R_u = %f\n', j_tmp, j_spr, R_u(j_spr,k_spr)); 
                                end 
                            end 
                        end
                        
                        % total_rest_length_free_edge = sum(R_v( min_idx:(j+min_idx), k_max)); 
                
                        % difference in points 
                        j_difference = next_comm_idx - (j + min_idx); 
                        
                        % then set that back with periodicity
                        j_reflected = mod(next_comm_idx + j_difference, j_max); 
                        
                        if j_reflected == 0
                            error('this shuold never be zero because zero is the comm point')
                        end 
                        
                        ring_point_reflected = X(:,j_reflected,1); 
                        
                        % leaflet rest height at reflected point
                        radial_stretch_j_reflected = radial_stretch_local(j_reflected, N_each, radial_stretch_center, extra_stretch_radial);
                        total_height_reflected = radial_stretch_j_reflected * sum(R_v(j_reflected, 1:(k_max-1)));
                        
                        

                    end 
                    
                    % relevant distances from each of three points 
                    % this is the intersection of three spheres 
                    F = @(p) [norm(ring_point - p) - total_height_current; norm(ring_point_reflected - p) - total_height_reflected; norm(comm_for_fused_edge - p) - total_rest_length_free_edge]; 
                                       
                    options = optimset('Display','off','TolFun',1e-20);
                    comm_interp_point = fsolve(F,[0;0;0],options);                    
                    
                    if (k==k_max) && debug_text 
                        fprintf("residual nonlinear solve for interp point %e = ", norm(F(comm_interp_point))); 
                    end 
                    
                    tangent = (comm_interp_point - ring_point); 
                    tangent = tangent / norm(tangent); 

                    % total radial rest length of this radial fiber 
                    total_rest_length = sum(R_v(j + min_idx, 1:(k-1))); 

                    % based on the rest length 
                    X(:,j + min_idx ,k) = total_rest_length * tangent * radial_stretch_center + ring_point;   
                    
                    if (k==k_max) && debug_text 
                        j
                        total_height_current
                        total_height_reflected
                        total_rest_length_free_edge
                        comm_interp_point
                        X_from_interpolation = X(:,j + min_idx ,k)
                        fprintf('\n')
                    end 
                    
                    if (j==4) 
                        'here'; 
                    end 
                    
                end 
            end 
        end 
    end 

    % error checking        
    tol = 1e-14; 
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

                if norm(X(:,j + min_idx,k_max) - X(:,j_reflected,k_max)) > tol 
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

function lambda_r = radial_stretch_local(j, N_each, radial_stretch_center, extra_stretch_radial)     

    points = [0 N_each/2 N_each]; 
    values = [extra_stretch_radial radial_stretch_center extra_stretch_radial]; 

    % lambda_r = interp1(points, values, mod(j,N_each)); 

    lambda_r = 1.0; radial_stretch_center; 
    
end 

