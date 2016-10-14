function [X] = build_reference_surface(leaflet)
%
% Builds reference coffee cone surface 
% 


N                       = leaflet.N; 
r                       = leaflet.filter.r; 
h                       = leaflet.filter.h; 
is_internal             = leaflet.is_internal; 
is_bc                   = leaflet.is_bc; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 


X           = zeros(3,j_max,k_max); 
X_flat      = nan * zeros(2,j_max,k_max); 

if (leaflet.min_angle < -pi) || (leaflet.max_angle > pi) 
    error('Outside allowable range of angles for current parameters'); 
end 


debug = false; 


if leaflet.radial_and_circumferential
   
    if isfield(leaflet, 'trapezoidal_flat_points')
        trapezoidal_flat_points = leaflet.trapezoidal_flat_points; 
    else 
        trapezoidal_flat_points = 0; 
    end 
    
    mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 
    ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

    % set the valve ring
    k = k_max; 
    for j=1:j_max
        
        X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), leaflet); 
        X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 

        % adjust final height so ring is in z = 0 plane 
        X(3,j,k) = X(3,j,k) - h; 

    end

    
    % Compute the "point" of the leaflet, which is the 1,1 index
    % This is where digonal fibers would connect if they were included 
    point = compute_point_of_leaflet(X_flat, leaflet); 
    
    spacing = norm( X_flat(:,j_max/2,k_max) -X_flat(:,j_max/2 + 1,k_max)); 
    
    % Running free edge point 
    % Leftmost is half spacing left from the point 
    % Plus spacing/2 for each extra leaflet point 
    current_x_flat = point + (trapezoidal_flat_points + 1)*[-spacing/2; 0];
    
    
    % Set free edges
   
    % first one is half an increment from the point 
    j = k_max; 
    k = 1;   
    X_flat(:,j,k) = current_x_flat; 
    
    % In 2d, work on a segment from the leftmost ring point to the 
    increment_left = (X_flat(:,1,k_max) - X_flat(:,j,k)) / (k_max - 1); 
    
    j = k_max - 1; 
    for k=2:(k_max - 1)
        % each one later gets the full increment 
        X_flat(:,j,k) = X_flat(:,j+1,k-1) + increment_left; 
        j = j - 1; 
    end 
    
    
    % flat part of free edge, if applicable 
    k = 1; 
    for j = (k_max+1):(k_max + trapezoidal_flat_points)
        X_flat(:,j,k) = X_flat(:,j-1,k) + [spacing; 0]; 
    end 
    
    
        
    j = k_max + 1 + trapezoidal_flat_points; 
    k = 1; 
    
    X_flat(:,j,k) = point + (trapezoidal_flat_points + 1)*[spacing/2; 0];
    
    increment_right = (X_flat(:,j_max,k_max) - X_flat(:,j,k)) / (k_max - 1); 
    
    j = k_max + 1 + trapezoidal_flat_points + 1; 
    for k=2:(k_max - 1)
        X_flat(:,j,k) = X_flat(:,j-1,k-1) + increment_right; 
        j = j+1; 
    end 

    if debug 
        figure; 
        plot(X_flat(1,:), X_flat(2,:), '-o')
        title('free edge and reference')
    end 
    
    
    
    % loop from free edge then up in k 
    j = k_max; 
    points_on_fiber = k_max; 
    for k=1:(k_max - 1)
        
        % nummber of points on current vertical fiber 
        % which decreases by one each index away from center of the leaflet 
        increment = (X_flat(:,j,k_max) - X_flat(:,j,k)) / (points_on_fiber - 1); 
        
        for k_tmp=(k+1):(k_max-1)
            X_flat(:,j,k_tmp) = X_flat(:,j,k_tmp - 1) + increment; 
        end
        
        j = j - 1;
        points_on_fiber = points_on_fiber - 1; 
    end 

    % from flat part up 
    k = 1; 
    for j = (k_max+1):(k_max + trapezoidal_flat_points)
        
        points_on_fiber = k_max; 
        increment = (X_flat(:,j,k_max) - X_flat(:,j,k)) / (points_on_fiber - 1); 
        
        for k_tmp=2:(k_max - 1)
            X_flat(:,j,k_tmp) = X_flat(:,j,k_tmp - 1) + increment; 
        end 
        
    end
    
    
    % from right free edge up 
    j = k_max + 1 + trapezoidal_flat_points ; 
    points_on_fiber = k_max; 
    
    for k=1:(k_max - 1)

        increment = (X_flat(:,j,k_max) - X_flat(:,j,k)) / (points_on_fiber - 1); 
        
        for k_tmp=(k+1):(k_max-1)
            X_flat(:,j,k_tmp) = X_flat(:,j,k_tmp - 1) + increment; 
        end
        
        j = j + 1;
        points_on_fiber = points_on_fiber - 1; 
    end
    
    
    if debug
        figure; 
        plot(X_flat(1,:), X_flat(2,:), '-o')
        title('full flat layout')
    end 
    
    
    % fill in 3d points 
    for k=1:k_max
        for j=1:j_max 
            if is_internal(j,k)
                X(:,j,k) = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet);
                X(3,j,k) = X(3,j,k) - h; 
            end 
        end 
    end 
    
else 
    
    mesh = linspace(leaflet.min_angle, leaflet.max_angle, N+1); 
    ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

    % set the valve ring 
    for j=1:j_max
        k = (N+2) - j; % ring coodinates 
        X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), leaflet); 
        X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 

        % adjust final height so ring is in z = 0 plane 
        X(3,j,k) = X(3,j,k) - h; 

    end 

    % fill in the 3d array 
    for j=1:j_max
        for k=1:k_max
            if leaflet.is_internal(j,k)
                X_flat(:,j,k)    = compute_intersection_diag_fibers(X_flat, j, k, leaflet); 
                X(:,j,k)         = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 
                X(3,j,k)         = X(3,j,k) - h; 
            end
        end 
    end

end 

