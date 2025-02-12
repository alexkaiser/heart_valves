function leaflets = import_valve_from_csv(valve, params, file_name, directory, run_inv_transform)

if exist('directory', 'var')
    file_with_path = fullfile(directory, file_name); 
else 
    file_with_path = file_name; 
end 

if ~exist('run_inv_transform', 'var')
    run_inv_transform = false; 
end 

[n_layers, n_leaflets] = size(params.layer_indices);

% just an allocation 
leaflets(n_layers, n_leaflets) = valve.leaflets(1);


for leaflet_num = 1:n_leaflets
    for layer_num = 1:n_layers 

        leaflet_temp = valve.leaflets(leaflet_num); 
        
        indices_global = params.layer_indices(layer_num, leaflet_num).indices_global;

        % file_with_path

        vertices = csvread(file_with_path)'; 

        if run_inv_transform
            vertices = coordinate_transformation_vertices(vertices, ...
                                                          valve.transformation_vertex_file, ...
                                                          valve.initial_rotation_aortic, ...
                                                          valve.initial_translation_aortic, ...
                                                          run_inv_transform); 
        end 

        j_max  = leaflet_temp.j_max; 
        k_max  = leaflet_temp.k_max; 

        for j = 1:j_max
            for k = 1:k_max
                vertex_idx = indices_global(j,k) + 1; 
                leaflet_temp.X(:,j,k) = vertices(:,vertex_idx);
            end 
        end 

        leaflets(layer_num, leaflet_num) = leaflet_temp; 
        
    end 
end 
    
    
    
    