function leaflet = aortic_free_edge_to_dirichlet_bc(leaflet)

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 

R_v    = leaflet.R_v; 

k = k_max; 

X = leaflet.X; 

debug = true; 

is_bc = leaflet.is_bc; 
linear_idx_offset         = zeros(j_max, k_max); 
point_idx_with_bc         = zeros(j_max, k_max); 

for comm_idx = 1:3

    % point one internal of commissure to point that m
    % N_each is a power of two 
    min_idx = (comm_idx-1)*N_each;         

    dj_interp = 1/N_each; 

    prev_comm_idx = min_idx; 
    if prev_comm_idx == 0
        prev_comm_idx = j_max; 
    end 
    
    comm_prev = X(:,prev_comm_idx,k); 
    comm_next = X(:,min_idx + N_each,k); 
    
    for j=1:(N_each-1)
        
        comm_interp_point = (1 - j*dj_interp) * comm_prev ...
                                 + j*dj_interp  * comm_next; 

        ring_point = X(:,j + min_idx ,1); 
        
        tangent = (comm_interp_point - ring_point); 
        tangent = tangent / norm(tangent); 
                             
        % total radial rest length of this radial fiber 
        total_rest_length = sum(R_v(j + min_idx, :)); 
        
        X(:,j + min_idx ,k) = total_rest_length * tangent + ring_point; 
        
        if is_bc(j + min_idx ,k)
            error('trying to set a bc as new position'); 
        end 
        
        is_bc(j + min_idx ,k) = true; 
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






