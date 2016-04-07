function X = build_reference_surface(filter_params)
%
% Builds reference coffee cone surface 
% 
% Mesh has a triangular layout 
% j==1 and k==1 are assumed to be the free edge 
% j+k == N+2 is the valve ring 
% 
% 
% 

N = filter_params.N; 
r = filter_params.r; 
h = filter_params.h; 

if (filter_params.min_angle < -pi) || (filter_params.max_angle > pi) 
    error('outside allowable range of angles for current parameters'); 
end 

X      = zeros(3,N+1,N+1); 
X_flat = zeros(2,N+1,N+1); 

mesh = linspace(filter_params.min_angle,filter_params.max_angle,N+1); 
ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

% set the valve ring 
for j=1:N+1
    k = (N+2) - j; % ring coodinates 
    X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
    X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
    
    % adjust final height so ring is in z = 0 plane 
    X(3,j,k) = X(3,j,k) - h; 
end 


% fill in the 3d array 
for j=1:N
    for k=1:N
        
        % in the triangle? 
        if ((j+k) < (N+2))
            X_flat(:,j,k) = compute_intersection(X_flat, j, k, filter_params); 
            X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
            X(3,j,k) = X(3,j,k) - h; 
        end

    end 
end


