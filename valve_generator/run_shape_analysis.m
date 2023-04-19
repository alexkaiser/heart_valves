

data_dir = '/Users/alex/data_to_remove/aortic_13673972_384_9c6b2e5_true_bicuspid_nostab_circ_3pt9_rad_1pt7_new_initial_cond/field_data'

file_name = 'aortic_no_partition_3840945.csv'

load aortic_no_partition_384_final_data.mat

indices_global = params.layer_indices(1).indices_global; 

run_inv_transform = true; 

basic_test = false; 
if basic_test 
    
    valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name, indices_global, data_dir, run_inv_transform); 

    valve_plot(valve_with_reference)

    [Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ...
     orifice_area_free_edge, orifice_area_all_points, free_edge_length_total] = shape_analysis_aortic(valve_with_reference.leaflets(1))


end 


run_all = false;



if run_all 

    
    times = csvread(fullfile(data_dir, 'times.txt')); 
    
    n_times = length(times); 

    leaflet = valve_with_reference.leaflets(1); 
    if isfield(leaflet, 'N_leaflets')
        N_leaflets = leaflet.N_leaflets; 
    else 
        N_leaflets = 3; 
    end 
    
    
    Gh = zeros(N_leaflets, n_times); 
    Eh = zeros(N_leaflets, n_times); 
    Eh_from_minimum = zeros(N_leaflets, n_times); 
    min_height_center = zeros(N_leaflets, n_times); 
    free_edge_length = zeros(N_leaflets, n_times); 
    orifice_area_free_edge = zeros(n_times, 1); 
    orifice_area_all_points = zeros(n_times, 1); 
    free_edge_length_total = zeros(n_times, 1); 
    
    % run down in index to preallocate struct array on first assignment 
    for i = length(times):-1:1
    
        % time files are zero indexed 
        file_name = sprintf('aortic_no_partition_384%04d.csv', i-1); 

        valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name, indices_global, data_dir, run_inv_transform); 

        [Gh(:,i), Eh(:,i), Eh_from_minimum(:,i), min_height_center(:,i), free_edge_length(:,i), ... 
         orifice_area_free_edge(i), orifice_area_all_points(i), free_edge_length_total(i)] = ...
            shape_analysis_aortic(valve_with_reference.leaflets(1)); 
    end

    save all_frames_data.mat
    
end

load_data = true; 
if load_data
    load all_frames_data.mat
end 

diameter = valve_with_reference.r * 2; 

% nondimensional parameters 
free_edge_length_over_diameter = free_edge_length / diameter;  
Gh_over_radius = Gh / valve_with_reference.r; 
Eh_over_Gh = Eh ./ Gh; 
Eh_from_min_over_Gh = Eh_from_minimum ./ Gh; 

width = 2.0; 

fig = figure; 
hold on;

plot(times, Gh(1,:), 'LineWidth',width); 
plot(times, Gh(2,:), 'LineWidth',width); 

plot(times, Eh(1,:), 'LineWidth',width); 
plot(times, Eh(2,:), 'LineWidth',width); 

plot(times, Eh_from_minimum(1,:), 'LineWidth',width); 
plot(times, Eh_from_minimum(2,:), 'LineWidth',width); 

% plot(times, min_height_center(1,:), 'LineWidth',width); 
% plot(times, min_height_center(2,:), 'LineWidth',width); 

legend('Gh 1', 'Gh 2', 'Eh 1', 'Eh 2', 'Eh from min 1', 'Eh from min 2')

set(fig, 'Position', [100, 100, 1000, 600])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'height_parameters.eps')

fig = figure; 
hold on;

plot(times, free_edge_length(1,:), 'LineWidth',width)
plot(times, free_edge_length(2,:), 'LineWidth',width)
legend('Free edge len 1', 'Free edge len 2')

set(fig, 'Position', [100, 100, 1000, 600])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'free_edge_length.eps')

% plot(times, free_edge_length_over_diameter(1,:), 'LineWidth',width)
% plot(times, free_edge_length_over_diameter(2,:), 'LineWidth',width)
% plot(times, free_edge_length_total(:))



fig = figure; 
hold on;

plot(times, orifice_area_free_edge, 'LineWidth',width)
plot(times, orifice_area_all_points, 'LineWidth',width)
legend('Area free edge', 'Area all points')

set(fig, 'Position', [100, 100, 1000, 600])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'area.eps')

fig = figure; 
hold on 
plot(times, free_edge_length_over_diameter(1,:), 'LineWidth',width)
plot(times, free_edge_length_over_diameter(2,:), 'LineWidth',width)

plot(times, Gh_over_radius(1,:), 'LineWidth',width)
plot(times, Gh_over_radius(2,:), 'LineWidth',width)

plot(times, Eh_over_Gh(1,:), 'LineWidth',width)
plot(times, Eh_over_Gh(2,:), 'LineWidth',width)

plot(times, Eh_from_min_over_Gh(1,:), 'LineWidth',width)
plot(times, Eh_from_min_over_Gh(2,:), 'LineWidth',width)

legend('Free edge over d 1', 'Free edge over d 2', 'Gh over r 1', 'Gh over r 2', 'Eh over Gh 1', 'Eh over Gh 2', 'Eh from min over Gh 1', 'Eh from min over Gh 2')

set(fig, 'Position', [100, 100, 1000, 600])
set(fig,'PaperPositionMode','auto')
printfig(fig, 'nondimensional_parameters.eps')


