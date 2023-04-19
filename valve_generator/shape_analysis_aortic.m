function [Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ...
          orifice_area_free_edge, orifice_area_all_points, free_edge_length_total] = shape_analysis_aortic(leaflet)

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 
X      = leaflet.X; 


if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    N_leaflets = 3; 
end 

Eh = zeros(N_leaflets,1); 
Gh = zeros(N_leaflets,1); 
min_height_center = zeros(N_leaflets,1); 
free_edge_length = zeros(N_leaflets,1);



% height 
for comm_idx = 1:N_leaflets

    min_idx = (comm_idx-1)*N_each;         
   
    center_idx_j = min_idx + N_each/2; 

    j = center_idx_j; 
    
    for k=2:k_max
    
        j_nbr_tmp = j;
        k_nbr_tmp = k-1; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = X(:,j,k);
        X_nbr = X(:,j_nbr,k_nbr); 

        Gh(comm_idx) = Gh(comm_idx) + norm(X_temp - X_nbr);        
    end 
    
    Eh(comm_idx) = X(3,j,k_max); 
    
    min_height_center(comm_idx) = min(X(3,j,:));
    
end

% measured from minimum of leaflet rather than center 
Eh_from_minimum = Eh - min_height_center; 

% free edge lengths 
for comm_idx = 1:N_leaflets

    min_idx = (comm_idx-1)*N_each;         
    
    for j=(1 + min_idx):(N_each + min_idx)
        k = k_max; 

        j_nbr_tmp = j-1; 
        k_nbr_tmp = k; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = X(:,j,k);
        X_nbr = X(:,j_nbr,k_nbr); 

        free_edge_length(comm_idx) = free_edge_length(comm_idx) + norm(X_temp - X_nbr);        
        
    end
    
end

free_edge_length_total = sum(free_edge_length); 

orifice_area_free_edge = polyarea(X(1,:,k_max), X(2,:,k_max)); 

% orifice area from conv hull 
% multiply all vectors by one over their norm 
% take convex hull then reinvert 
% for interior area 
vector_norms_xy = zeros(j_max,k_max); 

x_component_norm_inverted = zeros(j_max,k_max); 
y_component_norm_inverted = zeros(j_max,k_max); 

for j=1:j_max
    for k=1:k_max 
        vector_norms_xy(j,k) = norm(X(1:2,j,k)); 
        x_component_norm_inverted(j,k) = X(1,j,k) / vector_norms_xy(j,k)^2; 
        y_component_norm_inverted(j,k) = X(2,j,k) / vector_norms_xy(j,k)^2; 
    end 
end 

x_component = X(1,:,:); 
y_component = X(2,:,:); 
x_component_vector = x_component(:); 
y_component_vector = y_component(:); 

x_component_norm_inverted_vector = x_component_norm_inverted(:); 
y_component_norm_inverted_vector = y_component_norm_inverted(:); 

shrink_factor = 1; 
conv_hull_idx = boundary(x_component_norm_inverted_vector, y_component_norm_inverted_vector, shrink_factor); 

% figure; 
% plot(x_component_norm_inverted_vector, y_component_norm_inverted_vector, 'b*')
% hold on 
% plot(x_component_norm_inverted_vector(conv_hull_idx), y_component_norm_inverted_vector(conv_hull_idx), 'r-')

x_component_interior_polygon = x_component_vector(conv_hull_idx); 
y_component_interior_polygon = y_component_vector(conv_hull_idx); 

% figure; 
% plot(x_component_vector, y_component_vector, 'b*')
% hold on 
% plot(x_component_interior_polygon, y_component_interior_polygon, 'r-')

orifice_area_all_points = polyarea(x_component_interior_polygon, y_component_interior_polygon); 











