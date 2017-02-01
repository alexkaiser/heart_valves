function [X] = build_initial_fibers_bead_slip(leaflet)
%
% Builds reference coffee cone surface 
% 

N                       = leaflet.N; 
r                       = leaflet.r; 
is_internal             = leaflet.is_internal; 
is_bc                   = leaflet.is_bc; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 

left_papillary          = leaflet.left_papillary; 
right_papillary         = leaflet.right_papillary; 
free_edge_idx_left      = leaflet.free_edge_idx_left; 
free_edge_idx_right     = leaflet.free_edge_idx_right; 


X = NaN * zeros(3,j_max,k_max); 

debug = false; 

if leaflet.radial_and_circumferential
   
    % set the valve ring
    mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 
    X(:,:,k_max) = [r*cos(mesh); r*sin(mesh); zeros(size(mesh))]; 
    
    % Set free edge according to interpolating surface
    
    ring_l = X(:,1    ,k_max); 
    ring_r = X(:,j_max,k_max); 
    
    ds = 1/(j_max - 1); 
    
    interpolating_surf = @(s,t) t*(s*ring_r + (1-s)*ring_l) + (1-t)*(s*right_papillary + (1-s)*left_papillary);   
    t_of_s = @(s) abs(s-1/2) + 1/2 - ds; 
    
    % use the free edge arrays for indexing 
    for i=1:size(free_edge_idx_left, 1)

        j = free_edge_idx_left(i,1); 
        k = free_edge_idx_left(i,2); 

        s = (j-1)*ds; 
        X(:,j,k) = interpolating_surf(s,t_of_s(s)); 

    end
    
    for i=1:size(free_edge_idx_right, 1)

        j = free_edge_idx_right(i,1); 
        k = free_edge_idx_right(i,2); 

        s = (j-1)*ds; 
        X(:,j,k) = interpolating_surf(s,t_of_s(s)); 

    end      
    
    
    if debug 
        figure; 
        plot3(X(1,:), X(2,:), X(3,:), '-o'); 
        xlabel('x'); 
        ylabel('y'); 
        title('free edge and reference')
    end 
    
    % if true, takes crude guess at closed leaflet with curvature 
    % otherwise takes a linear interpolant from the free edge to the ring 
    pinched_interpolant = false; 
    
    if pinched_interpolant 
    
        % one dimensional mesh in straight line from commissure to commissure 
        
        ring_l = X(:,1    ,k_max); 
        ring_r = X(:,j_max,k_max);
        
        line_comm_to_comm = zeros(3,k_max+1); 
        ds = 1 / (j_max - 1); 
        for m=0:j_max
            line_comm_to_comm(:,m+1) = (m*ds)*ring_r + (1 - m*ds)*ring_l; 
        end 
        
        
        % fill in fibers interpolating between free edge and ring on each side 
        for i=1:size(free_edge_idx_left, 1)

            j = free_edge_idx_left(i,1); 
            k = free_edge_idx_left(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 
            
            X_line_commissure = 0.5 * (line_comm_to_comm(:,j) + X_free);
            
            num_points_first_half = floor(num_points/2); 
            ds_first = 1 / (num_points_first_half + 1); 
            
            num_points_second_half = num_points - num_points_first_half; 
            ds_second = 1 / (num_points_second_half + 1); 
            
            for m=1:num_points_first_half
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds_first)*X_line_commissure + (1 - m*ds_first)*X_free; 
            end
                   
            for m=1:num_points_second_half
                k_tmp = k + m + num_points_first_half; 
                X(:,j,k_tmp) = (m*ds_second)*X_ring + (1 - m*ds_second)*X_line_commissure;                      
            end

        end

        for i=1:size(free_edge_idx_right, 1)

            j = free_edge_idx_right(i,1); 
            k = free_edge_idx_right(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 
            
            X_line_commissure = 0.5 * (line_comm_to_comm(:,j) + X_free);

            num_points_first_half = floor(num_points/2); 
            ds_first = 1 / (num_points_first_half + 1); 
            
            num_points_second_half = num_points - num_points_first_half; 
            ds_second = 1 / (num_points_second_half + 1); 
            
            for m=1:num_points_first_half
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds_first)*X_line_commissure + (1 - m*ds_first)*X_free; 
            end
                   
            for m=1:num_points_second_half
                k_tmp = k + m + num_points_first_half; 
                X(:,j,k_tmp) = (m*ds_second)*X_ring + (1 - m*ds_second)*X_line_commissure;                      
            end
            
        end  
        
    else
    
        % linear interpolant from free edge 
        
        % fill in fibers interpolating between free edge and ring on each side 
        for i=1:size(free_edge_idx_left, 1)

            j = free_edge_idx_left(i,1); 
            k = free_edge_idx_left(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 

            for m=1:num_points
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds)*X_ring + (1 - m*ds)*X_free; 
            end 

        end

        for i=1:size(free_edge_idx_right, 1)

            j = free_edge_idx_right(i,1); 
            k = free_edge_idx_right(i,2); 

            % number of points on this fiber 
            num_points = k_max - k - 1; 

            % parameter spacing 
            ds = 1 / (k_max - k); 

            X_free = X(:,j,k); 
            X_ring = X(:,j,k_max); 

            for m=1:num_points
                k_tmp = k + m; 
                X(:,j,k_tmp) = (m*ds)*X_ring + (1 - m*ds)*X_free; 
            end

        end 
    
    end 
    
    if debug 
        figure; 
        
        x_component = squeeze(X(1,:,:)); 
        y_component = squeeze(X(2,:,:)); 
        z_component = squeeze(X(3,:,:)); 

        width = 1.5; 
        surf(x_component, y_component, z_component, 'LineWidth',width);
        
        axis equal 
        axis auto 
        
        xlabel('x'); 
        ylabel('y'); 
        title('leaflet')
    end 
    

else 
    error('diag not implemented with bead slip'); 
end 

