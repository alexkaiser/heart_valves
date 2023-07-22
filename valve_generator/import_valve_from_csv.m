function leaflet = import_valve_from_csv(valve, file_name, indices_global, directory, run_inv_transform)

if exist('directory', 'var')
    file_with_path = fullfile(directory, file_name); 
else 
    file_with_path = file_name; 
end 

if ~exist('run_inv_transform', 'var')
    run_inv_transform = false; 
end 

leaflet = valve.leaflets(1); 

% file_with_path

vertices = csvread(file_with_path)'; 

if run_inv_transform
    vertices = coordinate_transformation_vertices(vertices, ...
                                                  valve.transformation_vertex_file, ...
                                                  valve.initial_rotation_aortic, ...
                                                  valve.initial_translation_aortic, ...
                                                  run_inv_transform); 
end 

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 

for j = 1:j_max
    for k = 1:k_max
        vertex_idx = indices_global(j,k) + 1; 
        leaflet.X(:,j,k) = vertices(:,vertex_idx);
    end 
end 

