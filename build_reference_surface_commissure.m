function X = build_reference_surface_commissure(filter_params, left)
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

if left 
    papillary = [0; -filter_params.a; 0]; 
else 
    papillary = [0;  filter_params.a; 0]; 
end 

if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end 

X      = zeros(3,N+2,(N+3)/2); 
X_flat = zeros(2,N+2,(N+3)/2); 

mesh = linspace(filter_params.min_angle,filter_params.max_angle,N+3); 
ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

% set the valve ring 
for j=1:((N+3)/2)
    k = j; 
    X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
    X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
end 

for j=((N+3)/2 - 1):(N+2)
    k = N + 3 - j; 
    X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
    X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
end 

% fill in the 3d array 
for j=1:N+2
    for k=1:((N+3)/2)
        if is_internal_commissure(j,k,N)
            X_flat(:,j,k) = compute_intersection(X_flat, j, k, filter_params); 
            X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
        end

    end 
end





















