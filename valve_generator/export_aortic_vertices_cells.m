function export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform, export_cells, export_mechanics, export_coaptation)


[filepath,name,ext] = fileparts(file_name); 

name_no_number = name(1:end-4); 

n_layers = 3; 

n_leaflets = length(valve_with_reference.leaflets);

for leaflet_num = n_leaflets:-1:1
    for layer = n_layers:-1:1
        leaflets = import_valve_from_csv(valve_with_reference, params, file_name, data_dir, run_inv_transform); 
    end 
end 

cell_indices = zeros(8,1); 

j_max  = leaflets(1).j_max; 
k_max  = leaflets(1).k_max; 

n_vertices_total = j_max * k_max * n_layers * n_leaflets; 
vertices = zeros(3,n_vertices_total); 

if ~exist('export_mechanics', 'var')
    export_mechanics = false; 
end 

if ~exist('export_coaptation', 'var')
    export_coaptation = false; 
end 


indices_global_export = zeros(n_leaflets, n_layers, j_max, k_max); 

global_idx_running = 0; 

if export_mechanics
    
    layer_idx_ventricular = 1; 
    
    sigma_circ = zeros(n_vertices_total,1);
    sigma_rad = zeros(n_vertices_total,1);
    stress_circ = zeros(n_vertices_total,1);
    stress_rad = zeros(n_vertices_total,1);
    
    sigma_circ_ventricular  = zeros(n_vertices_total,1);
    sigma_rad_ventricular   = zeros(n_vertices_total,1);
    stress_circ_ventricular = zeros(n_vertices_total,1);
    stress_rad_ventricular  = zeros(n_vertices_total,1);
        
    for leaflet_num = 1:n_leaflets
                
        for layer = 1:n_layers

            [sigma_circ_temp, sigma_rad_temp, ~, ~, ~, stress_circ_temp, stress_rad_temp, ~, ~] ...
                = estimate_tangent_modulus_aortic_with_reference(leaflets(layer, leaflet_num), valve_with_reference.normal_thickness);

            leaflets(layer, leaflet_num).sigma_circ  = sigma_circ_temp; 
            leaflets(layer, leaflet_num).sigma_rad   = sigma_rad_temp; 
            leaflets(layer, leaflet_num).stress_circ = stress_circ_temp; 
            leaflets(layer, leaflet_num).stress_rad  = stress_rad_temp; 
           
            % ventricular layer for all leaflets 
            if layer == layer_idx_ventricular                
                for layer_temp = 1:n_layers
                    leaflets(layer_temp, leaflet_num).sigma_circ_ventricular  = sigma_circ_temp; 
                    leaflets(layer_temp, leaflet_num).sigma_rad_ventricular   = sigma_rad_temp; 
                    leaflets(layer_temp, leaflet_num).stress_circ_ventricular = stress_circ_temp; 
                    leaflets(layer_temp, leaflet_num).stress_rad_ventricular  = stress_rad_temp;                     
                end                 
            end 
                
        end
                
    end 
            
end 

if export_coaptation
    
    coapt_threshold = 0.2; 
    plots = false;
    
    coapt_distances = zeros(n_vertices_total,1);
    
    [~, ~, ~, ~, ~, ~, min_distances] = coaptation_analysis_aortic(leaflets, coapt_threshold, plots); 
    
    for leaflet_num = 1:n_leaflets                
        for layer = 1:n_layers
            leaflets(layer, leaflet_num).min_distances = squeeze(min_distances(leaflet_num,:,:));                 
        end                
    end 
        
end 

        
        
for leaflet_num = 1:n_leaflets
    for layer = 1:n_layers
        
        for j = 1:j_max
            for k = 1:k_max

                vertices(:,global_idx_running + 1) = leaflets(layer, leaflet_num).X(:,j,k); 

                if export_mechanics 
                    sigma_circ(global_idx_running + 1)  = leaflets(layer, leaflet_num).sigma_circ(j,k);
                    sigma_rad(global_idx_running + 1)   = leaflets(layer, leaflet_num).sigma_rad(j,k);
                    stress_circ(global_idx_running + 1) = leaflets(layer, leaflet_num).stress_circ(j,k);
                    stress_rad(global_idx_running + 1)  = leaflets(layer, leaflet_num).stress_rad(j,k);
                    
                    sigma_circ_ventricular(global_idx_running + 1)  = leaflets(layer, leaflet_num).sigma_circ_ventricular(j,k);
                    sigma_rad_ventricular(global_idx_running + 1)   = leaflets(layer, leaflet_num).sigma_rad_ventricular(j,k);
                    stress_circ_ventricular(global_idx_running + 1) = leaflets(layer, leaflet_num).stress_circ_ventricular(j,k);
                    stress_rad_ventricular(global_idx_running + 1)  = leaflets(layer, leaflet_num).stress_rad_ventricular(j,k);
                                                           
                end 
                
                if export_mechanics
                    coapt_distances(global_idx_running + 1) = leaflets(layer, leaflet_num).min_distances(j,k);
                end 
                                
                indices_global_export(leaflet_num, layer, j, k) = global_idx_running; 

                global_idx_running = global_idx_running + 1; 
            end 
        end 
    end 
end 



% fprintf("vertices capacity = %d\n", j_max * k_max * n_layers)
% fprintf("max idx placed = %d\n", global_idx_running)

% write needed vertices only as csv 
dlmwrite(strcat(name, '_vertices.csv'), vertices', 'delimiter', ' ', 'precision', 15); 


if export_mechanics
    filename_mechanics = sprintf('%s_mechanics.mat', name);    
    save(filename_mechanics, 'sigma_circ', 'sigma_rad', 'stress_circ', 'stress_rad', ...
        'sigma_circ_ventricular', 'sigma_rad_ventricular', 'stress_circ_ventricular', 'stress_rad_ventricular');
end 

if export_coaptation
    filename_mechanics = sprintf('%s_mechanics.mat', name);    
    save(filename_mechanics, 'coapt_distances', '-append');
end 

if export_cells

    cell_file_name = strcat(name_no_number, '_cells.csv');

    cell_file = fopen(cell_file_name, 'w'); 

    for leaflet_num = 1:n_leaflets
        for layer = 1:(n_layers-1)
            for j = 1:(j_max-1)
                for k = 1:(k_max-1)

                        [valid j_nbr_plus_one k_nbr] = get_indices(leaflets(layer, leaflet_num), j, k, j+1, k);                 
                        if ~valid
                            error('Found invalid index pair')
                        end 

                        [valid j_nbr k_nbr_plus_one] = get_indices(leaflets(layer, leaflet_num), j, k, j, k+1);                 
                        if ~valid
                            error('Found invalid index pair')
                        end                 


                        cell_indices(1) = indices_global_export(leaflet_num, layer  , j             , k             );
                        cell_indices(2) = indices_global_export(leaflet_num, layer+1, j             , k             );
                        cell_indices(3) = indices_global_export(leaflet_num, layer+1, j_nbr_plus_one, k             );
                        cell_indices(4) = indices_global_export(leaflet_num, layer  , j_nbr_plus_one, k             );
                        cell_indices(5) = indices_global_export(leaflet_num, layer  , j             , k_nbr_plus_one);
                        cell_indices(6) = indices_global_export(leaflet_num, layer+1, j             , k_nbr_plus_one);
                        cell_indices(7) = indices_global_export(leaflet_num, layer+1, j_nbr_plus_one, k_nbr_plus_one);
                        cell_indices(8) = indices_global_export(leaflet_num, layer  , j_nbr_plus_one, k_nbr_plus_one);

                        fprintf(cell_file, '8  ');                 
                        for cell_idx = 1:8
                            fprintf(cell_file, '%d  ', cell_indices(cell_idx)); 
                        end 
                        fprintf(cell_file, '\n'); 

                end 
            end 
        end 
    end 
    
    fclose(cell_file); 
end 


















