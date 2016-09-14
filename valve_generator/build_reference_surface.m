function [X is_internal is_bc]= build_reference_surface(leaflet)
%
% Builds reference coffee cone surface 
% 
% Mesh has a triangular layout 
% j==1 and k==1 are assumed to be the free edge 
% j+k == N+2 is the valve ring 
% 
% 
% 

N = leaflet.N; 
r = leaflet.filter.r; 
h = leaflet.filter.h; 

if (leaflet.min_angle < -pi) || (leaflet.max_angle > pi) 
    error('Outside allowable range of angles for current parameters'); 
end 

if ~leaflet.radial_and_circumferential
    
    X           = zeros(3,N+1,N+1); 
    X_flat      = zeros(2,N+1,N+1); 
    is_internal = zeros(N+1, N+1); 
    is_bc       = zeros(N+1, N+1); 

    mesh = linspace(leaflet.min_angle, leaflet.max_angle, N+1); 
    ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

    % set the valve ring 
    for j=1:N+1
        k = (N+2) - j; % ring coodinates 
        X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), leaflet); 
        X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 

        % adjust final height so ring is in z = 0 plane 
        X(3,j,k) = X(3,j,k) - h; 

        % ring points are boundary conditions 
        is_bc(j,k) = true; 
    end 

    % fill in the 3d array 
    for j=1:N
        for k=1:N
            % in the triangle? 
            if ((j+k) < (N+2))
                X_flat(:,j,k)    = compute_intersection(X_flat, j, k, leaflet); 
                X(:,j,k)         = cone_filter(X_flat(1,j,k), X_flat(2,j,k), leaflet); 
                X(3,j,k)         = X(3,j,k) - h; 
                is_internal(j,k) = true;  
            end
        end 
    end

else 
    error('Radial and circumferential fibers not implemented '); 
end 

