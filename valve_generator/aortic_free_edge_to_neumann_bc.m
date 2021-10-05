function leaflet = aortic_free_edge_to_neumann_bc(leaflet)

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 

is_bc = leaflet.is_bc; 

% turn all off 
is_bc(:,k_max) = false; 

linear_idx_offset         = zeros(j_max, k_max); 
point_idx_with_bc         = zeros(j_max, k_max); 

k=k_max; 
for j=(N_each * (1:3))
    is_bc(j,k) = true; 
end 

% pinches commissure points 
if isfield(leaflet, 'pinch_commissure') && leaflet.pinch_commissure
    if ~isfield(leaflet, 'N_to_pinch')
        error('must supply N_to_pinch if leaflet.pinch_commissure is true')
    end 
    
    for leaflet_idx=1:3

        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (leaflet_idx-1)*N_each;         

        prev_comm_idx = min_idx; 
        if prev_comm_idx == 0
            prev_comm_idx = j_max; 
        end 
        
        for j=1:leaflet.N_to_pinch
            j_current = j + min_idx; 
            is_bc(j_current, k_max) = true; 
            
            % zero indexed prev_comm_idx minus j, number past the comm 
            j_reflected_temp = mod(prev_comm_idx,j_max) - j; 
            % then set that back with periodicity
            j_reflected = mod(j_reflected_temp,j_max); 
            if j_reflected == 0
                error('this shuold never be zero because zero is the comm point')
            end 
            
            is_bc(j_reflected,k_max) = true; 
            
        end 
            
    end     
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
