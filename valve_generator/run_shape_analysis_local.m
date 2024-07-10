diary shape_analysis_results.txt

times_file = fopen('times.txt','r');
times = fscanf(times_file,'%f',[1 inf]); 

systole_time = 1.48416600000000; 
diastole_time = 1.28805000000000; 

[~, sys_idx_matlab] = min(abs(times - systole_time))
[~, dia_idx_matlab] = min(abs(times - diastole_time))

% frames are zero indexed 
sys_idx = sys_idx_matlab + 1 
dia_idx = dia_idx_matlab + 1 

n_times = 2; 
idx_frames = [775, 893]; 

time_printed = false; 

cycle_duration = 0.8; 
cycle_number = 2; 
t_min_pressure = 1.42; 
t_max_pressure = 1.52; 
t_min_flow = 1.25; 
t_mid_flow = 1.5; 

run_inv_transform = true; 


load('aortic_no_partition_384_final_data.mat', 'valve_with_reference', 'params');  

indices_global_1 = params.layer_indices(1).indices_global; 
indices_global_2 = params.layer_indices(2).indices_global; 
indices_global_3 = params.layer_indices(3).indices_global; 

% diastolic only 
% time files are zero indexed 
file_name_diastole = sprintf('aortic_no_partition_384%04d.csv', dia_idx); 

% current dir is data directory 
data_dir = cd; 

valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name_diastole, indices_global_1, data_dir, run_inv_transform); 

% layer copies 
leaflet_layer_2 = import_valve_from_csv(valve_with_reference, file_name_diastole, indices_global_2, data_dir, run_inv_transform);                 
leaflet_layer_3 = import_valve_from_csv(valve_with_reference, file_name_diastole, indices_global_3, data_dir, run_inv_transform); 


[Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ... 
    ~, ~, free_edge_length_total] = ...
    shape_analysis_aortic(valve_with_reference.leaflets(1), leaflet_layer_2, leaflet_layer_3); 


diameter = valve_with_reference.r * 2; 

% nondimensional parameters 
free_edge_length_over_diameter = free_edge_length / diameter;
Gh_over_radius = Gh / valve_with_reference.r;
Eh_over_Gh = Eh ./ Gh;
Eh_from_min_over_Gh = Eh_from_minimum ./ Gh;

coapt_plots = false; 
coapt_threshold = 0.2; 
[coapt_height, coapt_reserve_height, belly_height, coapt_min, coapt_max, fig] = coaptation_analysis_aortic(leaflet_layer_3, coapt_threshold, coapt_plots); 

r_vbr = leaflet_layer_3.skeleton.r; 
r_ic  = leaflet_layer_3.skeleton.r; 
h     = leaflet_layer_3.skeleton.normal_height;   
x = h - coapt_max; 
y = coapt_min; 

% eval_aortic_formulas(r_vbr, r_ic, h, belly_height, coapt_reserve_height, Gh(1), Eh(1), free_edge_length(1)/2, x, y); 

% fprintf('%s & %s & %s & ', names_struct(data_idx).circ_over_d, names_struct(data_idx).circ, names_struct(data_idx).rad); 
% if table_1
%     fprintf('%s & ', names_struct(data_idx).circ_over_d);
% elseif table_2 
%     fprintf('%s & ', names_struct(data_idx).comment); 
% else 
%     error('must pick table type')
% end 

fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & ', Gh(1), Gh_over_radius(1), min_height_center(1), Eh(1), Eh_from_minimum(1), ...
                                Eh_over_Gh(1), Eh_from_min_over_Gh(1), free_edge_length(1), free_edge_length_over_diameter(1), coapt_reserve_height); 

fprintf('\\\\ \n'); 
fprintf('\\hline \n'); 


    

file_name_systole = sprintf('aortic_no_partition_384%04d.csv', sys_idx); 

valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name_systole, indices_global_1, data_dir, run_inv_transform); 

% layer copies 
leaflet_layer_2 = import_valve_from_csv(valve_with_reference, file_name_systole, indices_global_2, data_dir, run_inv_transform);                 
leaflet_layer_3 = import_valve_from_csv(valve_with_reference, file_name_systole, indices_global_3, data_dir, run_inv_transform); 

[~, ~, ~, ~, ~, orifice_area_free_edge, orifice_area_all_points] = ...
    shape_analysis_aortic(valve_with_reference.leaflets(1), leaflet_layer_2, leaflet_layer_3); 


% flow related quantities 
[delta_p_mean, delta_p_max, q_full_cycle, q_systole, q_diff_full_systole, dt] = compute_gradient_stroke_vol(cycle_duration, cycle_number, t_min_pressure, t_max_pressure, t_min_flow, t_mid_flow); 

diameter = valve_with_reference.r * 2; 

% nondimensional parameters 
free_edge_length_over_diameter = free_edge_length / diameter;
Gh_over_radius = Gh / valve_with_reference.r;
Eh_over_Gh = Eh ./ Gh;
Eh_from_min_over_Gh = Eh_from_minimum ./ Gh;


% fprintf('%s & %s & %s & ', names_struct(data_idx).circ_over_d, names_struct(data_idx).circ, names_struct(data_idx).rad); 
% if table_1
%     fprintf('%s & ', names_struct(data_idx).circ_over_d);
% elseif table_2 
%     fprintf('%s & ', names_struct(data_idx).comment); 
% else 
%     error('must pick table type')
% end 

fprintf('%.2f & %.2f & %.2f  ', orifice_area_all_points, delta_p_mean, q_systole); 

fprintf('\\\\ \n'); 
fprintf('\\hline \n'); 

diary off 
