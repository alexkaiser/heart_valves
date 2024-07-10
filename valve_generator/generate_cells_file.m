
% add two default locations 
addpath ~/valve_generator 
addpath ~/mitral_fully_discrete/valve_generator 

data_dir = pwd; 
run_inv_transform = false; 

% if isfile('aortic_no_partition_384_final_data.mat')
%     load('aortic_no_partition_384_final_data.mat', 'valve_with_reference', 'params');  
% elseif isfile('aortic_no_partition_192_final_data.mat')
%     load('aortic_no_partition_192_final_data.mat', 'valve_with_reference', 'params');  
% end 

mat_file_list = dir('aortic_*final_data.mat'); 

if length(mat_file_list) ~= 1 
    error('found too many mat files')
end 

mat_file_name = mat_file_list(1).name;
load(mat_file_name , 'valve_with_reference', 'params');  

cell_file_exported = false; 

file_list = dir('aortic_*.csv'); 

for i = 1:length(file_list)
    
    file_name = file_list(i).name; 
    
    if ~(contains(file_name, '_vertices.csv') || contains(file_name, '_cells.csv') || contains(file_name, 'cylinder')) 

        if ~cell_file_exported
            export_cells = true;
        else
            export_cells = false;
        end
        
        file_name 
        export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform, export_cells); 
        
        if export_cells
            cell_file_exported = true; 
        end 

    end 
end 


% export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform); 




