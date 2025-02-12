function [Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ...
          orifice_area_free_edge, orifice_area_all_points, free_edge_length_total] = shape_analysis_aortic(leaflets, debug_plots)

j_max  = leaflets(1,1).j_max; 
k_max  = leaflets(1,1).k_max; 
N_each = leaflets(1,1).N_each; 

[n_layers, n_leaflets] = size(leaflets);

% all on first leaflet for now, no arrays
% Eh = zeros(N_leaflets,1); 
% Gh = zeros(N_leaflets,1); 
% min_height_center = zeros(N_leaflets,1); 
% free_edge_length = zeros(N_leaflets,1);

Eh = 0;
Gh = 0;
min_height_center = 0;
free_edge_length = 0;

layer_idx_ventricular = 1; 
layer_idx_aortic = n_layers; 


N_leaflets = 1; 

% height 

for leaflet_idx = 1:N_leaflets

    center_idx_j = N_each/2; 

    j = center_idx_j; 
    
    for k=2:k_max
    
        j_nbr_tmp = j;
        k_nbr_tmp = k-1; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflets(layer_idx_aortic, leaflet_idx), j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = leaflets(layer_idx_aortic, leaflet_idx).X(:,j,k);
        X_nbr = leaflets(layer_idx_aortic, leaflet_idx).X(:,j_nbr,k_nbr); 

        Gh(leaflet_idx) = Gh(leaflet_idx) + norm(X_temp - X_nbr);        
    end 
    
    Eh(leaflet_idx) = leaflets(layer_idx_aortic, leaflet_idx).X(3,j,k_max); 
    
    % min_height_center(comm_idx) = min(X(3,j,:));
    
    min_height_center(leaflet_idx) = min(min(leaflets(layer_idx_ventricular, leaflet_idx).X(3,:,:)));
    
end

% measured from minimum of leaflet rather than center 
Eh_from_minimum = Eh - min_height_center; 

% free edge lengths 
for leaflet_idx=1:N_leaflets
    
    for j=2:j_max
        k = k_max; 

        j_nbr_tmp = j-1; 
        k_nbr_tmp = k; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflets(layer_idx_ventricular, leaflet_idx), j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = leaflets(layer_idx_ventricular, leaflet_idx).X(:,j,k);
        X_nbr = leaflets(layer_idx_ventricular, leaflet_idx).X(:,j_nbr,k_nbr); 

        free_edge_length(leaflet_idx) = free_edge_length(leaflet_idx) + norm(X_temp - X_nbr);        
        
    end
    
end

free_edge_length_total = sum(free_edge_length); 


X_free_edge_points = [];
for leaflet_idx = 1:N_leaflets
    X_free_edge_points = [X_free_edge_points, leaflets(layer_idx_ventricular, leaflet_idx).X(:,:,k_max)];
end 

orifice_area_free_edge = polyarea(X_free_edge_points(1,:), X_free_edge_points(2,:)); 


x_component_vector = []; 
y_component_vector = []; 
x_component_norm_inverted_vector = []; 
y_component_norm_inverted_vector = []; 



for leaflet_idx = 1:n_leaflets
    for layer_idx = 1:n_layers 

        X_temp = leaflets(layer_idx, leaflet_idx).X; 
        
        % orifice area from conv hull 
        % multiply all vectors by one over their norm 
        % take convex hull then reinvert 
        % for interior area 
        vector_norms_xy = zeros(j_max,k_max); 

        x_component_norm_inverted = zeros(j_max,k_max); 
        y_component_norm_inverted = zeros(j_max,k_max); 

        for j=1:j_max
            for k=1:k_max 
                vector_norms_xy(j,k) = norm(X_temp(1:2,j,k)); 
                x_component_norm_inverted(j,k) = X_temp(1,j,k) / vector_norms_xy(j,k)^2; 
                y_component_norm_inverted(j,k) = X_temp(2,j,k) / vector_norms_xy(j,k)^2; 
            end 
        end 

        x_component = X_temp(1,:,:); 
        y_component = X_temp(2,:,:); 
        x_component_vector = [x_component_vector; x_component(:)]; 
        y_component_vector = [y_component_vector; y_component(:)]; 

        x_component_norm_inverted_vector = [x_component_norm_inverted_vector; x_component_norm_inverted(:)]; 
        y_component_norm_inverted_vector = [y_component_norm_inverted_vector; y_component_norm_inverted(:)]; 

    end 
end 
    
shrink_factor = 1; 
conv_hull_idx = boundary(x_component_norm_inverted_vector, y_component_norm_inverted_vector, shrink_factor); 

if exist('debug_plots', 'var') && debug_plots
    figure; 
    plot(x_component_norm_inverted_vector, y_component_norm_inverted_vector, 'b*')
    hold on           
    plot(x_component_norm_inverted_vector(conv_hull_idx), y_component_norm_inverted_vector(conv_hull_idx), 'r-')
end 

x_component_interior_polygon = x_component_vector(conv_hull_idx); 
y_component_interior_polygon = y_component_vector(conv_hull_idx); 

if exist('debug_plots', 'var') && debug_plots
    figure; 
    plot(x_component_vector, y_component_vector, 'b*')
    hold on 
    plot(x_component_interior_polygon, y_component_interior_polygon, 'r-')
end 

orifice_area_all_points = polyarea(x_component_interior_polygon, y_component_interior_polygon); 







