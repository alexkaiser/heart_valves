

run_inv_transform = true; 



basic_test = false; 
if basic_test 
    
    data_dir = '/Users/alex/data_to_remove/aortic_13673972_384_9c6b2e5_true_bicuspid_nostab_circ_3pt9_rad_1pt7_new_initial_cond/field_data'

    file_name = 'aortic_no_partition_3840945.csv'

    load aortic_no_partition_384_final_data.mat valve_with_reference params
    
    indices_global = params.layer_indices(1).indices_global; 
    
    close all 
    
    valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name, indices_global, data_dir, run_inv_transform); 

    valve_plot(valve_with_reference)

    [Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ...
     orifice_area_free_edge, orifice_area_all_points, free_edge_length_total] = shape_analysis_aortic(valve_with_reference.leaflets(1))


end 


two_frame_summary = true; 

if two_frame_summary
    
%     addpath('/Users/alex/Dropbox/NYU/research/mitral_fully_discrete/valve_generator')
    
%     dir_list = [... 
%         "aortic_14234920_384_a75f53c_true_bicuspid_circ_2pt8_1pt1d_rad_1pt4_new_initial_cond", ...
%         "aortic_14103984_384_28eb02c_true_bicuspid_circ_3pt0_1pt2d_rad_1pt7_new_initial_cond", ...
%         "aortic_14097127_384_4b0882c_true_bicuspid_circ_3pt5_1pt4d_rad_1pt7_new_initial_cond_2", ...
%         "aortic_14095958_384_67bdcde_true_bicuspid_circ_3pt9_1pt65d_rad_1pt7_new_initial_cond_2", ...    
%         "aortic_14149904_384_a75f53c_true_bicuspid_circ_4pt5_1pt8d_rad_1pt7_new_initial_cond", ...
%         "aortic_14728686_384_a75f53c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt7_circ_prestretch_1pt08", ...
%         "aortic_15051277_384_98d4f6c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt4", ...
%         "aortic_16698035_384_23a9213_true_bicuspid_circ_3pt9_1pt6d_rad_1pt65_less_bowl"]; 
   
    names_struct(1).dir =  "aortic_22425104_384_40405e1_true_bicuspid_circ_2pt7_1pt1d_rad_1pt6_semifinal_1"; 
    names_struct(1).circ =  "2.8"; 
    names_struct(1).circ_over_d =  "1.1d"; 
    names_struct(1).rad =  "1.4"; 
    names_struct(1).comment = ""; 

    names_struct(2).dir =  "aortic_22389012_384_38d07aa_true_bicuspid_circ_3pt0_1pt2d_rad_1pt7_semifinal_1"; 
    names_struct(2).circ =  "3.0"; 
    names_struct(2).circ_over_d =  "1.2d"; 
    names_struct(2).rad =  "1.7"; 
    names_struct(2).comment = "";
    
    names_struct(3).dir =  "aortic_22388339_384_c68c6b7__true_bicuspid_circ_3pt5_1pt4d_rad_1pt7_semifinal_1"; 
    names_struct(3).circ =  "3.5"; 
    names_struct(3).circ_over_d =  "1.4d"; 
    names_struct(3).rad =  "1.7"; 
    names_struct(3).comment = "";
    
    names_struct(4).dir =  "aortic_22387241_384_8e666fd_true_bicuspid_circ_3pt9_1pt65d_rad_1pt7_semifinal_1"; 
    names_struct(4).circ =  "3.9"; 
    names_struct(4).circ_over_d =  "1.6d"; 
    names_struct(4).rad =  "1.7"; 
    names_struct(4).comment = "At circ";
    
    names_struct(5).dir =  "aortic_22389594_384_8bb4d38_true_bicuspid_circ_4pt5_1pt8d_rad_1pt7_semifinal_1"; 
    names_struct(5).circ =  "4.5"; 
    names_struct(5).circ_over_d =  "1.8d"; 
    names_struct(5).rad =  "1.7"; 
    names_struct(5).comment = "";
    
%     names_struct(6).dir =  "aortic_14728686_384_a75f53c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt7_circ_prestretch_1pt08"; 
%     names_struct(6).circ =  "3.9"; 
%     names_struct(6).circ_over_d =  "1.57d"; 
%     names_struct(6).rad =  "1.7"; 
%     names_struct(6).comment = "Circ prestretch 1.08";
%     
%     names_struct(6).dir =  "aortic_15051277_384_98d4f6c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt4"; 
%     names_struct(6).circ =  "3.9"; 
%     names_struct(6).circ_over_d =  "1.57d"; 
%     names_struct(6).rad =  "1.4"; 
%     names_struct(6).comment = "Less tall";
%     
%     names_struct(7).dir =  "aortic_16698035_384_23a9213_true_bicuspid_circ_3pt9_1pt6d_rad_1pt65_less_bowl"; 
%     names_struct(7).circ =  "3.9"; 
%     names_struct(7).circ_over_d =  "1.6d"; 
%     names_struct(7).rad =  "1.65"; 
%     names_struct(7).comment = "Less bowl";
%     
%     names_struct(8).dir =  "aortic_21192946_384_ae92093_true_bicuspid_circ_3pt75_1pt5d_rad_1pt4_less_bowl_layers_2e4"; 
%     names_struct(8).circ =  "3.75"; 
%     names_struct(8).circ_over_d =  "1.5d"; 
%     names_struct(8).rad =  "1.4"; 
%     names_struct(8).comment = "Less bowl \& tall";
%     
%     names_struct(9).dir =  "aortic_21254567_384_fc9c4b3_true_bicuspid_circ_3pt0_1pt2d_rad_1pt4_less_bowl_layers_2e4"; 
%     names_struct(9).circ =  "3.0"; 
%     names_struct(9).circ_over_d =  "1.2d"; 
%     names_struct(9).rad =  "1.4"; 
%     names_struct(9).comment = "Less bowl \& tall";
    
    % data_idx_to_output = [4,6,7,8]; 
    data_idx_to_output = 1:5; 
    
    n_times = 2; 
    idx_frames = [775, 893]; 
    
    run_layers = false;
    if run_layers
        layer_idx_range = 1:3; 
    else 
        layer_idx_range = 1; 
    end 
        
    for layer = layer_idx_range
    
        
        
%         suffixes = [sprintf("_diastole_layer_%d_%d", layer, idx_frames(1)-1), ...
%                     sprintf("_systole__layer_%d_%d", layer, idx_frames(2)-1), ...
%                     sprintf("_diastole_layer_%d_%d", layer, idx_frames(3)-1), ...
%                     sprintf("_systole__layer_%d_%d", layer, idx_frames(4)-1)]; 

        for i = 1:2:length(idx_frames)
            suffixes = [sprintf("_diastole_layer_%d_%d", layer, idx_frames(i)  -1), ...
                        sprintf("_systole__layer_%d_%d", layer, idx_frames(i+1)-1)]; 
        end 

        for i = 1:length(idx_frames)

            time_printed = false; 

            for data_idx = data_idx_to_output % 1:length(names_struct)

                dir = names_struct(data_idx).dir;  

                cd(convertStringsToChars(dir)); 

                data_dir = pwd; 

                load(fullfile(data_dir, 'aortic_no_partition_384_final_data.mat'), 'valve_with_reference', 'params');  
                times = csvread(fullfile(data_dir, 'times.txt')); 

                indices_global = params.layer_indices(layer).indices_global; 
                
                % proceed if this timestep is not available 
%                 if idx_frames(i) > length(times)
%                     cd .. 
%                     continue; 
%                 end 
                
                % times_tmp = times(idx_frames)'; 

                if ~time_printed 
                    fprintf('time = %f\n', times(idx_frames(i)))
                    time_printed = true; 
                end 

                leaflet = valve_with_reference.leaflets(1); 
                if isfield(leaflet, 'N_leaflets')
                    N_leaflets = leaflet.N_leaflets; 
                else 
                    N_leaflets = 3; 
                end 

    %             Gh = zeros(N_leaflets); 
    %             Eh = zeros(N_leaflets); 
    %             Eh_from_minimum = zeros(N_leaflets); 
    %             min_height_center = zeros(N_leaflets); 
    %             free_edge_length = zeros(N_leaflets); 
    %             orifice_area_free_edge = zeros(n_times, 1); 
    %             orifice_area_all_points = zeros(n_times, 1); 
    %             free_edge_length_total = zeros(n_times, 1); 

                % time files are zero indexed 
                file_name = sprintf('aortic_no_partition_384%04d.csv', idx_frames(i)); 

                valve_with_reference.leaflets(1) = import_valve_from_csv(valve_with_reference, file_name, indices_global, data_dir, run_inv_transform); 

                [Gh, Eh, Eh_from_minimum, min_height_center, free_edge_length, ... 
                 orifice_area_free_edge, orifice_area_all_points, free_edge_length_total] = ...
                    shape_analysis_aortic(valve_with_reference.leaflets(1)); 

                run_stretch_plots = false; 
                if run_stretch_plots

                    fiber_output    = true; 
                    fiber_stride    = 8; 
                    stride_offset_j = 4; 

                    circ  = true; 
                    rad   = false; 

                    az = 0; 
                    el = 60; 

                    min_plot_cap = 1.0; 
                    max_plot_cap = 1.2; 

                    fig = figure; 
                    set(fig, 'Renderer', 'Painters');

                    [lambda_circ, lambda_rad, lambda_circ_mean, lambda_rad_mean, fig]  = compute_stretch_aortic_with_reference(valve_with_reference.leaflets(1), fig, fiber_stride, stride_offset_j, circ, rad, min_plot_cap, max_plot_cap);
                    view(az,el);

                    outname_circ = strcat('circ_stretch_aortic', suffixes(i));                 
                    print(fig, '-depsc', outname_circ);

                    fig = figure; 
                    set(fig, 'Renderer', 'Painters');
                    circ  = false; 
                    rad   = true; 
                    ratio = false; 
                    min_plot_cap = 1.0; 
                    max_plot_cap = 2.0; 
                    [lambda_circ, lambda_rad, lambda_circ_mean, lambda_rad_mean, fig]  = compute_stretch_aortic_with_reference(valve_with_reference.leaflets(1), fig, fiber_stride, stride_offset_j, circ, rad, min_plot_cap, max_plot_cap);
                    view(az,el);
                    outname_rad = strcat('rad_stretch_aortic', suffixes(i));                 
                    print(fig, '-depsc', outname_rad);

                    close all 
                    
                    % 'pause'
                end 



                diameter = valve_with_reference.r * 2; 

    %             Gh
    %             Eh
    %             Eh_from_minimum
    %             min_height_center
    %             free_edge_length
    %             orifice_area_free_edge 
    %             orifice_area_all_points 
    %             free_edge_length_total 

                % nondimensional parameters 
                free_edge_length_over_diameter = free_edge_length / diameter;
                Gh_over_radius = Gh / valve_with_reference.r;
                Eh_over_Gh = Eh ./ Gh;
                Eh_from_min_over_Gh = Eh_from_minimum ./ Gh;


                % fprintf('%s & %s & %s & ', names_struct(data_idx).circ_over_d, names_struct(data_idx).circ, names_struct(data_idx).rad); 
                fprintf('%s & %s & %s & ', names_struct(data_idx).circ_over_d); 

                fprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f ', orifice_area_all_points, Gh(1), Gh_over_radius(1), min_height_center(1), Eh(1), Eh_from_minimum(1), ...
                                                Eh_over_Gh(1), Eh_from_min_over_Gh(1), free_edge_length(1), free_edge_length_over_diameter(1)); 

                fprintf('\\\\ \n'); 
                fprintf('\\hline \n'); 


                cd .. 


            end


        end 
    
    end 
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

load_data = false; 
if load_data
    load all_frames_data.mat
end 


run_plots = false; 
if run_plots 
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

end 



run_stl_export = false; 
if run_stl_export
       
    names_struct(1).dir =  "aortic_14234920_384_a75f53c_true_bicuspid_circ_2pt8_1pt1d_rad_1pt4_new_initial_cond"; 
    names_struct(1).circ =  "2.8"; 
    names_struct(1).circ_over_d =  "1.1d"; 
    names_struct(1).rad =  "1.4"; 
    names_struct(1).comment = ""; 

    names_struct(2).dir =  "aortic_14103984_384_28eb02c_true_bicuspid_circ_3pt0_1pt2d_rad_1pt7_new_initial_cond"; 
    names_struct(2).circ =  "3.0"; 
    names_struct(2).circ_over_d =  "1.2d"; 
    names_struct(2).rad =  "1.7"; 
    names_struct(2).comment = "";
    
    names_struct(3).dir =  "aortic_14097127_384_4b0882c_true_bicuspid_circ_3pt5_1pt4d_rad_1pt7_new_initial_cond_2"; 
    names_struct(3).circ =  "3.5"; 
    names_struct(3).circ_over_d =  "1.4d"; 
    names_struct(3).rad =  "1.7"; 
    names_struct(3).comment = "";
    
    names_struct(4).dir =  "aortic_14095958_384_67bdcde_true_bicuspid_circ_3pt9_1pt65d_rad_1pt7_new_initial_cond_2"; 
    names_struct(4).circ =  "3.9"; 
    names_struct(4).circ_over_d =  "1.57d"; 
    names_struct(4).rad =  "1.7"; 
    names_struct(4).comment = "Circumference";
    
    names_struct(5).dir =  "aortic_14149904_384_a75f53c_true_bicuspid_circ_4pt5_1pt8d_rad_1pt7_new_initial_cond"; 
    names_struct(5).circ =  "4.5"; 
    names_struct(5).circ_over_d =  "1.8d"; 
    names_struct(5).rad =  "1.7"; 
    names_struct(5).comment = "";
    
%     names_struct(6).dir =  "aortic_14728686_384_a75f53c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt7_circ_prestretch_1pt08"; 
%     names_struct(6).circ =  "3.9"; 
%     names_struct(6).circ_over_d =  "1.57d"; 
%     names_struct(6).rad =  "1.7"; 
%     names_struct(6).comment = "Circ prestretch 1.08";
%     
    names_struct(6).dir =  "aortic_15051277_384_98d4f6c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt4"; 
    names_struct(6).circ =  "3.9"; 
    names_struct(6).circ_over_d =  "1.57d"; 
    names_struct(6).rad =  "1.4"; 
    names_struct(6).comment = "Less tall";
    
    names_struct(7).dir =  "aortic_16698035_384_23a9213_true_bicuspid_circ_3pt9_1pt6d_rad_1pt65_less_bowl"; 
    names_struct(7).circ =  "3.9"; 
    names_struct(7).circ_over_d =  "1.6d"; 
    names_struct(7).rad =  "1.65"; 
    names_struct(7).comment = "Less bowl";
    
    names_struct(8).dir =  "aortic_21192946_384_ae92093_true_bicuspid_circ_3pt75_1pt5d_rad_1pt4_less_bowl_layers_2e4"; 
    names_struct(8).circ =  "3.75"; 
    names_struct(8).circ_over_d =  "1.5d"; 
    names_struct(8).rad =  "1.4"; 
    names_struct(8).comment = "Less bowl \& tall";
    
    names_struct(9).dir =  "aortic_21254567_384_fc9c4b3_true_bicuspid_circ_3pt0_1pt2d_rad_1pt4_less_bowl_layers_2e4"; 
    names_struct(9).circ =  "3.0"; 
    names_struct(9).circ_over_d =  "1.2d"; 
    names_struct(9).rad =  "1.4"; 
    names_struct(9).comment = "Less bowl \& tall";
    
    
    
%     dir_list = ["aortic_14095958_384_67bdcde_true_bicuspid_circ_3pt9_1pt65d_rad_1pt7_new_initial_cond_2"]; 

    n_times = 2; 
    idx_frames = [242, 407]; 

    suffixes = [sprintf("_diastole_%d", idx_frames(1)-1), sprintf("_systole_%d", idx_frames(2)-1)]; 
    
    % run down in index to preallocate struct array on first assignment 
    for i = 1:length(idx_frames)
    
        time_printed = false; 
        
        for data_idx = 1:length(names_struct)

            dir = names_struct(data_idx).dir;  

            cd(convertStringsToChars(dir)); 

            data_dir = pwd; 

            load(fullfile(data_dir, 'aortic_no_partition_384_final_data.mat'), 'valve_with_reference', 'params');  
            times = csvread(fullfile(data_dir, 'times.txt')); 

            indices_global = params.layer_indices(1).indices_global; 
            

            times_tmp = times(idx_frames)'; 

            if ~time_printed 
                fprintf('time = %f\n', times_tmp(i))
                time_printed = true; 
            end 
            
            leaflet = valve_with_reference.leaflets(1); 
            if isfield(leaflet, 'N_leaflets')
                N_leaflets = leaflet.N_leaflets; 
            else 
                N_leaflets = 3; 
            end 

            % time files are zero indexed 
            file_name = sprintf('aortic_no_partition_384%04d.csv', idx_frames(i)-1); 
            
            export_aortic_vertices_cells(file_name, valve_with_reference, params, data_dir, run_inv_transform); 
            
            cd .. 
            
        end 
    end     
end 





