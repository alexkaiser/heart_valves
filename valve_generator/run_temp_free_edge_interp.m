function [free_edge_len, free_edge_interp_points] = run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max_from_center, power)

% runs the interpolation via input curve 
% and returns the resulting length 

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 

R_u    = leaflet.R_u; 
R_v    = leaflet.R_v; 

if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    error("must provide N_leaflets")
end 

% if N_leaflets ~= 2
%     error("N_leaflets must be equal to 2")
% end 

if ~exist('extra_stretch_radial', 'var')
    extra_stretch_radial = 1.0; 
end

if ~exist('power', 'var')
    power = 4.0; 
end

X = leaflet.X; 


debug_transform = false; 

coord_transform = true; 
if coord_transform
    
    tol = 1e-14; 
    
    comm_1 = X(:,1,k_max); 
    comm_2 = X(:,j_max,k_max); 
    
    % translate comm1 to origin in xy 
    translation_comm_1_to_origin = [-X(1,1,k_max); -X(2,1,k_max); 0];
    
    comm_1 = comm_1 + translation_comm_1_to_origin;
    comm_2 = comm_2 + translation_comm_1_to_origin;
    X = X + translation_comm_1_to_origin;
    
    % put comm_2  to the y axis     
    rotation = rotation_matrix_z(pi - atan2(comm_2(2), comm_2(1)));
    
    comm_1 = rotation * comm_1;
    comm_2 = rotation * comm_2;
    
    for k=1:k_max
        X(:,:,k) = rotation * X(:,:,k);
    end 
    
    if abs(comm_2(2)) > tol 
        error('comm_2 did not land on x axis');
    end 
    
    translation_center_on_y = [-comm_2(1)/2; 0; 0];
    comm_1 = comm_1 + translation_center_on_y;
    comm_2 = comm_2 + translation_center_on_y;
    X = X + translation_center_on_y;
    
    if abs(comm_1(1) + comm_2(1)) > tol
        error('x components not centered on y')
    end 
        
    if debug_transform
%         figure; 
%         x_component = squeeze(X(1,:,:)); 
%         y_component = squeeze(X(2,:,:)); 
%         z_component = squeeze(X(3,:,:)); 
% 
%         width = 1.0; 
%         surf(x_component, y_component, z_component, 'LineWidth',width);
% 
%         axis equal 
%         axis auto 
%         hold on 
    end 
    
end 



% start with initial free edge 
% commissure points stay fixed 
free_edge_interp_points = X(:,:,k_max); 


for j=2:N_each

    k = k_max; 

    ring_point = X(:,j,1); 

    th = atan2(ring_point(2), ring_point(1));  

    % make a litle 
    strained_len_total = extra_stretch_radial * sum(R_v(j, :)); 

    % cm apart at middle 
%         y_free_edge_end = y_max_from_center * sign(sin(th)) * sin(th)^2; 
    % power = 4; 
    y_free_edge_end = y_max_from_center * sign(sin(th)) * abs(sin(th)^power);
    % y_free_edge_end = y_max_from_center * ring_point(2);
    % y_free_edge_end = 0; 
    % this would put the two free edges exactly coinciding 

    % if using exact x 
    % then (y_diff^2 + height^2) = strained_len_total^2 
    % so height is given as 
    interp_height = sqrt(strained_len_total^2 - (ring_point(2) - y_free_edge_end)^2) ; 

    % comm_interp_point = [ring_point(1) ; (r/2) * sin(th); comm_prev(3)];                 
    free_edge_interp_points(:,j) = [ring_point(1) ; y_free_edge_end; interp_height + ring_point(3)]; 

    % total radial rest length of this radial fiber 
    total_rest_length = sum(R_v(j, 1:(k-1))); 

    tangent = (free_edge_interp_points(:,j) - ring_point); 
    tangent = tangent / norm(tangent); 

    % based on the rest length 
    X(:,j,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

end 

free_edge_len = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u); 



if coord_transform
    free_edge_interp_points = free_edge_interp_points - translation_center_on_y; 
    free_edge_interp_points = rotation' * free_edge_interp_points;
    free_edge_interp_points = free_edge_interp_points - translation_comm_1_to_origin;    
end


if debug_transform
    figure; 
    x_component = squeeze(X(1,:,:)); 
    y_component = squeeze(X(2,:,:)); 
    z_component = squeeze(X(3,:,:)); 
    hold on 

    width = 1.0; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

    x_component = squeeze(free_edge_interp_points(1,:,:)); 
    y_component = squeeze(free_edge_interp_points(2,:,:)); 
    z_component = squeeze(free_edge_interp_points(3,:,:)); 

    width = 1.0; 
    plot3(x_component, y_component, z_component, 'r-o', 'LineWidth',width);

    
    axis equal 
    axis auto 
    hold on 
end 


