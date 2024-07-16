
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
linear_idx_offset         = zeros(j_max, k_max); 
point_idx_with_bc         = zeros(j_max, k_max); 


if N_leaflets == 2 
    free_edge_interp_points = find_free_edge_interp_points_true_bicuspid(leaflet, extra_stretch_radial, extra_stretch_circ); 
end 



for comm_idx = 1:N_leaflets

    % point one internal of commissure to point that m
    % N_each is a power of two 
    min_idx = (comm_idx-1)*N_each;         

    dj_interp = 1/N_each; 

    prev_comm_idx = min_idx; 
    if prev_comm_idx == 0
        prev_comm_idx = j_max; 
    end 
    
    comm_prev = X(:,prev_comm_idx,k_max); 
    comm_next = X(:,min_idx + N_each,k_max); 
    
    for j=1:(N_each-1)
        for k=k_range
        
            ring_point = X(:,j + min_idx ,1); 
            
            % total radial rest length of this radial fiber 
            total_rest_length = sum(R_v(j + min_idx, 1:(k-1))); 
            
            if N_leaflets == 2 
                comm_interp_point = free_edge_interp_points(:,j + min_idx); 
            else 
                comm_interp_point = (1 - j*dj_interp) * comm_prev ...
                                       + j*dj_interp  * comm_next; 
            end 


            tangent = (comm_interp_point - ring_point); 
            tangent = tangent / norm(tangent); 


            % based on the rest length 
            X(:,j + min_idx ,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

            % based on putting the free edge as a triangle between commn points 
            % generally bad 
            % X(:,j + min_idx ,k) = (k/k_max) * (comm_interp_point - ring_point) + ring_point; 
            
            if is_bc(j + min_idx ,k)
                error('trying to set a bc as new position'); 
            end 

            if k == k_max
                is_bc(j + min_idx ,k) = true; 
            end 
        end 
    end 

end 

% % allows a bit of equalizing at the free edge by comms 
% % looks terrible 
% if (N_leaflets == 2)
%     for comm_idx = 1:N_leaflets
%         % point one internal of commissure to point that m
%         % N_each is a power of two 
%         neumann_fraction_each_half = 1/8;                 
%         min_idx = (comm_idx-1)*N_each;         
%         for j=1:(N_each-1)
%             if abs(j - (N_each/2)) > (N_each * (1/2 - neumann_fraction_each_half))
%                 is_bc(j + min_idx ,k_max) = false;  
%             end 
%         end 
%     end 
% end 

% 
% if (N_leaflets == 2)
%     for comm_idx = 1:N_leaflets
%         % point one internal of commissure to point that m
%         % N_each is a power of two 
%         neumann_fraction_center = 1/4;                 
%         min_idx = (comm_idx-1)*N_each;         
%         for j=1:(N_each-1)
%             if abs(j - (N_each/2)) < (N_each * neumann_fraction_center)
%                 is_bc(j + min_idx ,k_max) = false;  
%             end 
%         end 
%     end 
% end 

debug_lengths = true; 
if debug_lengths

%     free_edge_length_single_loaded = 0; 
%     free_edge_length_single_rest = 0; 
%     for j=1:N_each
%         k=k_max; 
% 
%         j_nbr_tmp = j-1; 
%         k_nbr_tmp = k; 
%         [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
%         if ~valid 
%             error('trying to compute lengths with an invalid rest length')
%         end
% 
%         X_temp = X(:,j,k);
%         X_nbr = X(:,j_nbr,k_nbr); 
% 
%         j
%         loaded_len = norm(X_temp - X_nbr)
%         rest_len = R_u(j_spr,k_spr)
%         
%         free_edge_length_single_loaded = free_edge_length_single_loaded + norm(X_temp - X_nbr);
%         
%         free_edge_length_single_rest = free_edge_length_single_rest + R_u(j_spr,k_spr); 
%     end     
%     
%     free_edge_length_single_loaded
%     free_edge_length_single_rest
    
    'after reinterpolating free edge'
    [free_edge_length_single_loaded, free_edge_length_single_rest] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
end 




leaflet.X = X; 

% for j=1:j_max 
%     if ~is_bc(j,k)
%         error('did not set all bcs')
%     end 
% end 

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






