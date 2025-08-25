function leaflet = aortic_free_edge_to_dirichlet_bc(leaflet, extra_stretch_radial, extra_stretch_circ)

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 

R_u    = leaflet.R_u; 
R_v    = leaflet.R_v; 

if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    N_leaflets = 3; 
end 
    


full_leaflet_interp = true; 
if full_leaflet_interp
    k_range = 2:k_max; 
else
    k_range = k_max; 
end 

if ~exist('extra_stretch_radial', 'var')
    extra_stretch_radial = 1.0; 
end

if ~exist('extra_stretch_circ', 'var')
    extra_stretch_circ = 1.0; 
end

if N_leaflets ~= 2
    if extra_stretch_circ ~= 1.0 
        error('extra circ stretch only implemented for true bicuspid')
    end 
end 


X = leaflet.X; 

debug = true; 

is_bc = leaflet.is_bc; 
linear_idx_offset = zeros(j_max, k_max); 
point_idx_with_bc = zeros(j_max, k_max); 


free_edge_equalization = true; 

if free_edge_equalization
    free_edge_interp_points = find_free_edge_interp_points_aortic(leaflet, extra_stretch_radial, extra_stretch_circ); 
else 
    dj_interp = 1/N_each; 

    % always boundaries 
    comm_prev = X(:,1,k_max); 
    comm_next = X(:,j_max,k_max); 
end 


for j=2:N_each
    for k=k_range

        ring_point = X(:,j,1); 

        % total radial rest length of this radial fiber 
        total_rest_length = sum(R_v(j, 1:(k-1))); 

        if free_edge_equalization
            comm_interp_point = free_edge_interp_points(:,j); 
        else 
            % old interpolation style on chord 
            comm_interp_point = (1 - (j-1)*dj_interp) * comm_prev ...
                                   + (j-1)*dj_interp  * comm_next; 

        end 

        tangent = (comm_interp_point - ring_point); 
        tangent = tangent / norm(tangent); 

        % based on the rest length 
        X(:,j,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

        if is_bc(j,k)
            error('trying to set a bc as new position'); 
        end 

        if k == k_max
            is_bc(j,k) = true; 
        end 
    end 
end 


debug_lengths = true; 
if debug_lengths
    'after reinterpolating free edge'
    [free_edge_length_single_loaded, free_edge_length_single_rest] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
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
    title('In aortic free edge to dirichlet bc')
end 

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






