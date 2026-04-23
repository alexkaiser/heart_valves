
% file_name = 'meshes/aortic_f258d7a_trileaflet_4_1_26/aortic_no_partition_384_final_data.mat'; 

file_name = 'aortic_no_partition_384_final_data.mat'; 

load(file_name, 'valve', 'valve_with_reference', 'N')

close all 

N_leaflets = valve.leaflets.N_leaflets;


% annulus figure 
annulus_figure = true; 
if annulus_figure 

    fig = figure; 
    set(fig, 'Renderer', 'Painters');

    x_temp = []; 
    y_temp = [];
    z_temp = [];

    j_max = valve_with_reference.leaflets(1).j_max;
    k_max = valve_with_reference.leaflets(1).k_max;

    for i = 1:N_leaflets

        x_temp = [x_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(1,1,:)))']; 
        y_temp = [y_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(2,1,:)))']; 
        z_temp = [z_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(3,1,:)))']; 

        x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,:,1))]; 
        y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,:,1))]; 
        z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,:,1))];         

        x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,j_max,:))']; 
        y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,j_max,:))']; 
        z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,j_max,:))']; 

    end 

    x_temp = [x_temp(end), x_temp]; 
    y_temp = [y_temp(end), y_temp]; 
    z_temp = [z_temp(end), z_temp];     
        
    width = 3; 
    plot3(x_temp, y_temp, z_temp, 'k', 'LineWidth', width)
    hold on 

    axis equal 
    axis off 

    az = -122.6192;
    el = 23.2580;
    view(az,el)

    camproj('perspective')
    
    set(fig, 'Position', [0, 0, 1000, 500])
    set(fig,'PaperPositionMode','auto')
    print(fig, '-depsc', 'aortic_annulus');

    % fig = valve_plot(valve_with_reference, fig);
    % 
    % print(fig, '-depsc', 'aortic_initial_cond');
    % 
end 




tangent_mod_plots = false; 
if tangent_mod_plots 
    
    % constitutive law version 
    valve_with_reference_no_solve = valve; 

    % kill off the old leaflet structure, new one has different fields, 
    % which makes matlab complain about assigning it to a structure array 
    valve_with_reference_no_solve  = rmfield(valve_with_reference_no_solve , 'leaflets'); 

    valve_with_reference_no_solve.leaflets(1) = set_rest_lengths_and_constants_aortic(valve.leaflets(1), valve); 
    
    fiber_output    = true; 
    fiber_stride    = 8; 
    stride_offset_j = 4; 

    circ  = true; 
    rad   = false; 
    ratio = false; 

    az = -90; 
    el =   0; 
        
    max_plot_cap = 4e8; 
        
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    
    [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, fig]  = estimate_tangent_modulus_aortic_with_reference(valve_with_reference_no_solve.leaflets(1), valve.normal_thickness, fig, fiber_stride, stride_offset_j, circ, rad, ratio, max_plot_cap);

    hold on 

    

    x_temp = []; 
    y_temp = [];
    z_temp = [];

    i = 1
    x_temp = [x_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(1,1,:)))']; 
    y_temp = [y_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(2,1,:)))']; 
    z_temp = [z_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(3,1,:)))']; 

    x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,:,1))]; 
    y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,:,1))]; 
    z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,:,1))];         

    x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,j_max,:))']; 
    y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,j_max,:))']; 
    z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,j_max,:))'];      
    

    

        
    width = 6; 
    plot3(x_temp, y_temp, z_temp, 'k', 'LineWidth', width)

    view(az,el);
    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic');


    az = -122.6192;
    el = 23.2580;
    view(az,el)

    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_annulus_view');


    az = -211.0944
    el = 10.7652

    view(az,el)
    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_side_view');

end 



tangent_mod_plots_reference = true; 
if tangent_mod_plots_reference
    
    fiber_output    = true; 
    fiber_stride    = 8; 
    stride_offset_j = 4; 

    circ  = true; 
    rad   = false; 
    ratio = false; 

    az = -90; 
    el =   0; 
        
    max_plot_cap = 4e8; 
        
    fig = figure; 
    set(fig, 'Renderer', 'Painters');

    x_temp = []; 
    y_temp = [];
    z_temp = [];

    j_max = valve_with_reference.leaflets(1).j_max;
    k_max = valve_with_reference.leaflets(1).k_max;

    
    for i=1:N_leaflets
        [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, fig]  = estimate_tangent_modulus_aortic_with_reference(valve_with_reference.leaflets(i), valve.normal_thickness, fig, fiber_stride, stride_offset_j, circ, rad, ratio, max_plot_cap);
    
        hold on 
    
        x_temp = [x_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(1,1,:)))']; 
        y_temp = [y_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(2,1,:)))']; 
        z_temp = [z_temp, fliplr(squeeze(valve_with_reference.leaflets(i).X(3,1,:)))']; 
    
        x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,:,1))]; 
        y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,:,1))]; 
        z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,:,1))];         
    
        x_temp = [x_temp, squeeze(valve_with_reference.leaflets(i).X(1,j_max,:))']; 
        y_temp = [y_temp, squeeze(valve_with_reference.leaflets(i).X(2,j_max,:))']; 
        z_temp = [z_temp, squeeze(valve_with_reference.leaflets(i).X(3,j_max,:))'];      
            
       

    end 

    width = 6; 
    plot3(x_temp, y_temp, z_temp, 'k', 'LineWidth', width)

    view(az,el);
    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_tri');


    az = -122.6192;
    el = 23.2580;
    view(az,el)

    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_annulus_view_tri');


    az = -211.0944
    el = 10.7652

    view(az,el)
    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_side_view_tri');

    az = -180;
    el = 90;

    view(az,el)
    print(fig, '-depsc', 'circ_tangent_mod_aortic_schematic_top_tri');


end 


% ic_plot = false; 
% if ic_plot
% 
%     valve_plot(valve_with_reference);
% 
% 
%     axis equal;
% 
% end