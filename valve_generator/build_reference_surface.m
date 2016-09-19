function [X] = build_reference_surface(leaflet)
%
% Builds reference coffee cone surface 
% 
% Mesh has a triangular layout 
% j==1 and k==1 are assumed to be the free edge 
% j+k == N+2 is the valve ring 
% 
% 

N           = leaflet.N; 
r           = leaflet.filter.r; 
h           = leaflet.filter.h; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc; 
j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 

X           = zeros(3,j_max,k_max); 
X_flat      = nan * zeros(2,j_max,k_max); 

if (leaflet.min_angle < -pi) || (leaflet.max_angle > pi) 
    error('Outside allowable range of angles for current parameters'); 
end 

if leaflet.radial_and_circumferential
   
    mesh = linspace(leaflet.min_angle, leaflet.max_angle, N); 
    ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

    % set the valve ring
    k = k_max; 
    for j=1:j_max
        
        X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), leaflet); 
        X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 

        % adjust final height so ring is in z = 0 plane 
        X(3,j,k) = X(3,j,k) - h; 

    end
    
    figure; 
    plot(X_flat(1,:), X_flat(2,:), '-o'); 
    title('ring only')
    
    
    % Compute the "point" of the leaflet, which is the 1,1 index
    % This is where digonal fibers would connect if they were included 
    point = compute_point_of_leaflet(X_flat, leaflet); 
    
    spacing = norm( X_flat(:,j_max/2,k_max) -X_flat(:,j_max/2 + 1,k_max)); 
    
    
    % Set free edges
   
    % first one is half an increment from the point 
    j = j_max/2; 
    k = 1;   
    X_flat(:,j,k) = point + [-spacing/2; 0]; 
    
    % In 2d, work on a segment from the leftmost ring point to the 
    increment_left = (X_flat(:,1,k_max) - X_flat(:,j,k)) / (k_max - 1); 
    
    j = j_max/2 - 1; 
    for k=2:(k_max - 1)
        % each one later gets the full increment 
        X_flat(:,j,k) = X_flat(:,j+1,k-1) + increment_left; 
        j = j - 1; 
    end 
    
        
    j = j_max/2 + 1; 
    k = 1; 
    
    X_flat(:,j,k) = point + [spacing/2; 0];
    
    increment_right = (X_flat(:,j_max,k_max) - X_flat(:,j,k)) / (k_max - 1); 
    
    j = N/2 + 2; 
    for k=2:(k_max - 1)
        X_flat(:,j,k) = X_flat(:,j-1,k-1) + increment_right; 
        j = j+1; 
    end 
    
    
    % loop from free edge then up in k 
    j = j_max/2; 
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

    j = j_max/2 + 1; 
    points_on_fiber = k_max; 
    
    for k=1:(k_max - 1)

        increment = (X_flat(:,j,k_max) - X_flat(:,j,k)) / (points_on_fiber - 1); 
        
        for k_tmp=(k+1):(k_max-1)
            X_flat(:,j,k_tmp) = X_flat(:,j,k_tmp - 1) + increment; 
        end
        
        j = j + 1;
        points_on_fiber = points_on_fiber - 1; 
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
                X_flat(:,j,k)    = compute_intersection(X_flat, j, k, leaflet); 
                X(:,j,k)         = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 
                X(3,j,k)         = X(3,j,k) - h; 
            end
        end 
    end

end 

