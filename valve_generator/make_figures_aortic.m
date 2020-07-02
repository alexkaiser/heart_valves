
file_name = 'meshes/aortic_semifinal_6fa5a93/aortic_384_final_data.mat'; 

load(file_name, 'valve', 'valve_with_reference', 'N')

sigma_circ_mean = valve_with_reference.leaflets(1).sigma_circ_mean 
sigma_rad_mean  = valve_with_reference.leaflets(1).sigma_rad_mean 

mean_ratio_average_first = sigma_circ_mean / sigma_rad_mean

fprintf('mean_ratio_average_first = %.14f', mean_ratio_average_first)

sigma_circ_clark = 5.8e7; 

rel_err_clark = abs((sigma_circ_mean - sigma_circ_clark)/sigma_circ_clark) 


closed_valve_figure = true; 
if closed_valve_figure 
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    fiber_stride    = 8; 
    stride_offset_j = 4; 
    exploded_view = false;
    fig = surf_plot_fiber_subset(valve.leaflets(1), fig, fiber_stride, stride_offset_j, exploded_view); 
    
    axis equal 
    axis tight 
    axis off 
    
    print(fig, '-depsc', 'aortic_closed_diag_view');
    
    view(0,90)
    print(fig, '-depsc', 'aortic_closed_top_view');
    
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    fiber_stride    = 8; 
    stride_offset_j = 4; 
    exploded_view = false;
    one_leaflet = true; 
    fig = surf_plot_fiber_subset(valve.leaflets(1), fig, fiber_stride, stride_offset_j, exploded_view, one_leaflet); 
    axis equal 
    axis tight 
    axis off 
    print(fig, '-depsc', 'aortic_closed_one_leaflet_diag_view');
end 

closed_valve_exploded = false; 
if closed_valve_exploded
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    fiber_stride    = 8; 
    stride_offset_j = 4; 
    exploded_view = true; 
    fig = surf_plot_fiber_subset(valve.leaflets(1), fig, fiber_stride, stride_offset_j, exploded_view); 
    view(0,90); 
    print(fig, '-depsc', 'aortic_closed_top_view_exploded'); 
end 



% annulus figure 
annulus_figure = false; 
if annulus_figure 

    fig = figure; 
    set(fig, 'Renderer', 'Painters');
        
    for subplot_num = 1:2

        subplot(1,2,subplot_num)
        x_temp = squeeze(valve.leaflets(1).X(1,:,1)); 
        y_temp = squeeze(valve.leaflets(1).X(2,:,1)); 
        z_temp = squeeze(valve.leaflets(1).X(3,:,1)); 
        x_temp = [x_temp(end), x_temp]; 
        y_temp = [y_temp(end), y_temp]; 
        z_temp = [z_temp(end), z_temp]; 

        width = 3; 

        plot3(x_temp, y_temp, z_temp, 'k', 'LineWidth', width)
        hold on 

        for j_temp = (valve.leaflets(1).N_each * (1:3))

            x_temp = squeeze(valve.leaflets(1).X(1,j_temp,:)); 
            y_temp = squeeze(valve.leaflets(1).X(2,j_temp,:)); 
            z_temp = squeeze(valve.leaflets(1).X(3,j_temp,:)); 
            plot3(x_temp, y_temp, z_temp, 'k', 'LineWidth', width)

        end 

        axis equal 
        axis off 
        camproj('perspective')
        % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])

        if subplot_num == 1
            view(2,0)
        else 
            view(-4,28)
        end 
        
    end 
    
    set(fig, 'Position', [0, 0, 1000, 500])
    set(fig,'PaperPositionMode','auto')
    print(fig, '-depsc', 'aortic_annulus_two_panel');
    
end 

tension_figure = false; 
if tension_figure

    fiber_output    = true; 
    fiber_stride    = 8; 
    stride_offset_j = 4; 

    height_plot = false; 

    circ  = true; 
    rad   = false; 
    ratio = false; 

    az = -90; 
    el =   0; 
        
    leaflet = valve.leaflets(1); 

    x_size = 500; 
    y_size = 500; 
    
    fig = figure; 
    set(fig, 'Renderer', 'Painters');

    colorbar_on = false; 
    
    total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig, colorbar_on)
    view(az,el);
    print(fig, '-depsc', 'circ_tension_aortic');
    
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    circ  = false; 
    rad   = true; 
    ratio = false; 
    colorbar_on = false; 
    total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig, colorbar_on)
    view(az,el);
    print(fig, '-depsc', 'rad_tension_aortic');

    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    circ  = false; 
    rad   = false; 
    ratio = true; 
    colorbar_on = false; 
    total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig, colorbar_on)
    view(az,el);
    print(fig, '-depsc', 'ratio_tension_aortic');
    
    % plots tension plots, don't look great without modification... 
    % [fig_circ fig_rad fig_ratio] = make_aortic_plots(valve.leaflets(1)); 

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
    view(az,el);
    print(fig, '-depsc', 'circ_tangent_mod_aortic');
    
    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    circ  = false; 
    rad   = true; 
    ratio = false; 
    [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, fig]  = estimate_tangent_modulus_aortic_with_reference(valve_with_reference_no_solve.leaflets(1), valve.normal_thickness, fig, fiber_stride, stride_offset_j, circ, rad, ratio, max_plot_cap);
    view(az,el);
    print(fig, '-depsc', 'rad_tangent_mod_aortic');

    fig = figure; 
    set(fig, 'Renderer', 'Painters');
    circ  = false; 
    rad   = false; 
    ratio = true; 
    max_plot_cap_ratio = 60; 
    [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, fig]  = estimate_tangent_modulus_aortic_with_reference(valve_with_reference_no_solve.leaflets(1), valve.normal_thickness, fig, fiber_stride, stride_offset_j, circ, rad, ratio, max_plot_cap_ratio);
    view(az,el);
    print(fig, '-depsc', 'ratio_tangent_mod_aortic');

end 
        















