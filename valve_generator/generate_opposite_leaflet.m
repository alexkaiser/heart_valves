function leaflet = generate_opposite_leaflet(leaflet_current)
%
% Generates a leaflet which is opposite the current leaflet  
% and shares a free edge 
% 
% Free edge points are set to be NOT internal 
% 

leaflet.N           = leaflet_current.N; 
leaflet.total_angle = 2*pi - leaflet_current.total_angle; 
leaflet.min_angle   = pi + leaflet.total_angle/2.0;
leaflet.max_angle   = pi - leaflet.total_angle/2.0; 

leaflet.r = leaflet_current.r; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = leaflet_current.radial_and_circumferential; 

if ~leaflet.radial_and_circumferential
    error('generate opposite only implemented for radial and circumferential')
end 

leaflet.j_max               = leaflet_current.j_max; 
leaflet.k_max               = leaflet_current.k_max; 
leaflet.free_edge_idx_left  = leaflet_current.free_edge_idx_left;  
leaflet.free_edge_idx_right = leaflet_current.free_edge_idx_right; 
leaflet.chordae_idx_left    = leaflet_current.chordae_idx_left;  
leaflet.chordae_idx_right   = leaflet_current.chordae_idx_right; 
leaflet.is_internal         = leaflet_current.is_internal; 
leaflet.is_bc               = leaflet_current.is_bc; 
leaflet.linear_idx_offset   = leaflet_current.linear_idx_offset; 
leaflet.point_idx_with_bc   = leaflet_current.point_idx_with_bc; 

% Spring constants in two directions 
leaflet.alpha    = leaflet_current.alpha; 
leaflet.beta     = leaflet_current.beta; 
leaflet.p_0      = leaflet_current.p_0;
leaflet.ref_frac = leaflet_current.ref_frac; 

% No chordae here 
leaflet.chordae_tree = false; 

% Pull from other leaflet 
X = NaN * zeros(size(leaflet_current.X)); 

% Unpack, pull free edges from other leaflet 
r                       = leaflet.r; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 
free_edge_idx_left      = leaflet.free_edge_idx_left; 
free_edge_idx_right     = leaflet.free_edge_idx_right; 

% set the valve ring
mesh = linspace(leaflet.min_angle, leaflet.max_angle, j_max); 
X(:,:,k_max) = [r*cos(mesh); r*sin(mesh); zeros(size(mesh))]; 


for i=1:size(free_edge_idx_left, 1)
    j = free_edge_idx_left(i,1); 
    k = free_edge_idx_left(i,2); 
    X(:,j,k) = leaflet_current.X(:,j,k); 
    
    % free edges NOT internal on reflected leaflet 
    leaflet.is_internal(j,k) = false; 
end

for i=1:size(free_edge_idx_right, 1)
    j = free_edge_idx_right(i,1); 
    k = free_edge_idx_right(i,2); 
    X(:,j,k) = leaflet_current.X(:,j,k); 
    
    % free edges NOT internal on reflected leaflet 
    leaflet.is_internal(j,k) = false; 
end 


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

leaflet.X = X; 







