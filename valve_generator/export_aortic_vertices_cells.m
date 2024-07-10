function export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform, export_cells)


[filepath,name,ext] = fileparts(file_name); 

name_no_number = name(1:end-4); 

n_layers = 3; 

for layer = n_layers:-1:1
    leaflets(layer) = import_valve_from_csv(valve_with_reference, file_name, params.layer_indices(layer).indices_global, data_dir, run_inv_transform); 
end 

cell_indices = zeros(8,1); 

j_max  = leaflets(1).j_max; 
k_max  = leaflets(1).k_max; 

n_vertices_total = j_max * k_max * n_layers; 
vertices = zeros(3,n_vertices_total); 

indices_global_export = zeros(n_layers, j_max, k_max); 

global_idx_running = 0; 

for layer = 1:n_layers
    for j = 1:j_max
        for k = 1:k_max

            vertices(:,global_idx_running + 1) = leaflets(layer).X(:,j,k); 

            indices_global_export(layer,j,k) = global_idx_running; 
            
            global_idx_running = global_idx_running + 1; 
        end 
    end 
end 

% fprintf("vertices capacity = %d\n", j_max * k_max * n_layers)
% fprintf("max idx placed = %d\n", global_idx_running)

% write needed vertices only as csv 
dlmwrite(strcat(name, '_vertices.csv'), vertices', 'delimiter', ' ', 'precision', 15); 

if export_cells

    cell_file_name = strcat(name_no_number, '_cells.csv');

    cell_file = fopen(cell_file_name, 'w'); 

    for layer = 1:(n_layers-1)
        for j = 1:j_max
            for k = 1:(k_max-1)

                    [valid j_nbr_plus_one k_nbr] = get_indices(leaflets(layer), j, k, j+1, k);                 
                    if ~valid
                        error('Found invalid index pair')
                    end 

                    [valid j_nbr k_nbr_plus_one] = get_indices(leaflets(layer), j, k, j, k+1);                 
                    if ~valid
                        error('Found invalid index pair')
                    end                 


                    cell_indices(1) = indices_global_export(layer  , j             , k  );
                    cell_indices(2) = indices_global_export(layer+1, j             , k  );
                    cell_indices(3) = indices_global_export(layer+1, j_nbr_plus_one, k  );
                    cell_indices(4) = indices_global_export(layer  , j_nbr_plus_one, k  );
                    cell_indices(5) = indices_global_export(layer  , j             , k_nbr_plus_one);
                    cell_indices(6) = indices_global_export(layer+1, j             , k_nbr_plus_one);
                    cell_indices(7) = indices_global_export(layer+1, j_nbr_plus_one, k_nbr_plus_one);
                    cell_indices(8) = indices_global_export(layer  , j_nbr_plus_one, k_nbr_plus_one);

                    fprintf(cell_file, '8  ');                 
                    for cell_idx = 1:8
                        fprintf(cell_file, '%d  ', cell_indices(cell_idx)); 
                    end 
                    fprintf(cell_file, '\n'); 

            end 
        end 
    end 

    fclose(cell_file); 
end 


















