
% add default location
addpath ~/heart_valves/valve_generator 

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
        
        [filepath,name,ext] = fileparts(file_name); 

        name_no_number = name(1:end-4); 

        frame_num = str2num(name(end-3:end));
        
        if (frame_num == 662) || (frame_num == 837) 
            export_mechanics = true 
            export_coaptation = true
        else 
            export_mechanics = false
            export_coaptation = false
        end 
        
        export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform, export_cells, export_mechanics, export_coaptation); 
        
        if export_cells
            cell_file_exported = true; 
        end 

    end 
end 


% export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform); 




