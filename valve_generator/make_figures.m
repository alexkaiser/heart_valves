% 
% script to output many figures for thesis 
% 


schematics = false; 

if schematics 

    N = 128; 

    radial       = true; 
    bead_slip    = true; 
    attached     = false; 
    leaflet_only = false; 
    optimization = false; 
    decreasing_tension = true; 

    valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, decreasing_tension); 

%    output_leaflet_mesh_schematic(valve); 
%
    output_mesh_schematic(valve); 
%
%     output_valve_ring(valve); 
% 
%     output_dec_tension_curve(); 
% 
     N_jacobian = 32; 
%     output_jacobian_figure(N_jacobian); 
    
end 



patch_schematic = false; 
if patch_schematic
    load /Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/mitral_tree_128_final_data
    
    j_range_anterior = valve.leaflets(1).j_range_anterior; 
    
    X_patch = valve.leaflets(1).X(:,(32-10:32+10), 48:60); 
    
    x_component = squeeze(X_patch(1,:,:)); 
    y_component = squeeze(X_patch(2,:,:)); 
    z_component = squeeze(X_patch(3,:,:)); 
    width = 1.0; 
    
    fig = figure; 
    surf(x_component, y_component, z_component, 'LineWidth',width);
        
    axis equal 
    axis auto 
    grid off 
    axis off 
    view(-150,16)
    printfig(fig, 'leaflet_patch')
    
end     
    


static_plots = false; 
if static_plots  

    load /Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/mitral_tree_256_final_data

    fig = figure; 
    valve_plot(valve,fig)
    printfig(fig, 'valve_diagonal')

    view(0,90)
    printfig(fig, 'valve_top')

    view(0,0)
    printfig(fig, 'valve_side')

    view(90,50)
    printfig(fig, 'valve_mouth')
    
    view(90,0)
    printfig(fig, 'valve_front')

    'stop'

    grid off 
    view(-118,20)
    axis([-2.5 0.0 .3 .7 -2.1 -1.9]);
    zoom(2)
    axis off 
    printfig(fig, 'tree_detail')
    close all; 

    fig = figure; 
    valve_plot(valve_with_reference,fig)
    printfig(fig, 'valve_with_ref_diagonal')

    view(0,90)
    printfig(fig, 'valve_with_ref_top')

    view(0,0)
    printfig(fig, 'valve_with_ref_side')

    view(90,50)
    printfig(fig, 'valve_with_ref_mouth')
   
    
    grid off 
    view(90,0)
    axis equal
    axis([-2.5 0.0 -2 2 -2.8 -1.6]);
    zoom(2)
    axis off 
    printfig(fig, 'valve_with_ref_front')
end 


static_plots_comm = false; 
    
if static_plots_comm   

    load /Users/alex/Dropbox/NYU/research/mitral_fully_discrete/valve_generator/meshes/plot_meshes/7_25_17_comm_with_80mmhg_new_scaling_512_e265671/mitral_tree_256_final_data

    fig = figure; 
    valve_plot(valve,fig)
    printfig(fig, 'valve_diagonal_comm')

    view(0,90)
    printfig(fig, 'valve_top_comm')

    view(0,0)
    printfig(fig, 'valve_side_comm')

    view(90,50)
    printfig(fig, 'valve_mouth_comm')
    
    view(90,0)
    printfig(fig, 'valve_front_comm')
    
    'stop'
    
    grid off 
    view(-118,20)
    axis([-2.5 0.0 .3 .7 -2.1 -1.9]);
    zoom(2)
    axis off 
    printfig(fig, 'tree_detail_comm')
    close all; 
    
    fig = figure; 
    valve_plot(valve_with_reference,fig)
    printfig(fig, 'valve_with_ref_diagonal_comm')
    
    view(0,90)
    printfig(fig, 'valve_with_ref_top_comm')
    
    view(0,0)
    printfig(fig, 'valve_with_ref_side_comm')

    view(90,50)
    printfig(fig, 'valve_with_ref_mouth_comm')

    
    grid off 
    view(90,0)
    axis equal
    axis([-1.8 0.0 -1.8 1.8 -2.8 -1]);
    zoom(2)
    axis([-2.3 0.0 -1.4 1.4 -2.7 -1.18]);
    axis off 
    printfig(fig, 'valve_with_ref_front_comm')
end 



tension_plots = false; 
if tension_plots

    debug = true; 
    if debug
        n = 32; 
    else
        n = 512;
    end 
    
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
        
    circ = false; 
    radial = true; 
    fig_anterior_rad = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    % savefig(fig_anterior_rad, 'anterior_tension_plot_radial')
    % close(fig_anterior)
    view(90,0)
    axis equal; 
    % axis off; 
    anterior_axes_radial = gca
    limits_horiz_view = axis; 
    
    
    
    
    circ = true; 
    radial = false; 
    fig_anterior_circ = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    % savefig(fig_anterior_circ, 'anterior_tension_plot_circ')
    view(90,0)
    axis equal; 
    % axis off; 
    anterior_axes_circ = gca 
    limits_horiz_tmp = axis; 
    limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 
    % close(fig_anterior)
    
    
%     circ = true; 
%     radial = true; 
%     fig_anterior_both = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
%     % savefig(fig_anterior_both, 'anterior_tension_plot_both')
%     % close(fig_anterior)
%     axis equal; 
%     axis off; 
%     anterior_axes_both = gca; 

    
%     grid off 
%     view(-118,20)
%     axis([-2.5 0.0 .3 .7 -2.1 -1.9]);
%     zoom(2)
%     axis off 
%     printfig(fig_anterior_both, 'anterior_tension_tree_detail')
%     
%     grid off 
%     view(90,0)
%     axis equal
%     axis([-2.5 0.0 -2 2 -2.8 -1.6]);
%     zoom(2)
%     axis off 
%     printfig(fig_anterior_both, 'anterior_tension_free_edge_detail')
    
    anterior = false; 
    circ = true; 
    radial = false; 
    fig_posterior_circ = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    % savefig(fig_posterior, 'posterior_tension_plot_circ')
    % close(fig_posterior)
    view(90,0)
    axis equal;
    axis off; 
    posterior_axes_circ = gca 
    limits_horiz_tmp = axis; 
    limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

    circ = false; 
    radial = true; 
    fig_posterior_rad = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    % savefig(fig_posterior, 'posterior_tension_plot_radial')
    % close(fig_posterior)
    view(90,0)
    axis equal;
    axis off; 
    posterior_axes_rad = gca
    limits_horiz_tmp = axis; 
    limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 
    

%     circ = true; 
%     radial = true; 
%     fig_posterior_both = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
%     % savefig(fig_posterior, 'posterior_tension_plot_both')
%     % close(fig_posterior)
%     axis equal;
%     axis off; 
%     posterior_axes_both = gca; 
%     limits_horiz_tmp = axis; 
%     limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

    % figure out the min and max of each axis
    % y_mins = anterior_axes_circ.Ymin;
    % set to be largest window on all four head-on plots
    % also need to make sure that all are evenly zoomed 
    
    figure(fig_anterior_circ);

    % figure out and manually set the figure output size in pixels 
    lims = axis; 
    x_size = lims(2)-lims(1); 
    y_size = lims(4)-lims(3); 
    x_pixels = floor(100 * x_size); 
    y_pixels = floor(100 * y_size);
    
    
    axis(limits_horiz_view); 
    set(fig_anterior_circ, 'Position', [100, 100, x_pixels, y_pixels])
    set(fig_anterior_circ, 'PaperPosition', [100, 100, x_pixels, y_pixels])
    %set(fig_anterior_circ,'PaperPositionMode','auto')
    axis off; 
    % printfig(fig_anterior_circ, 'anterior_tension_plot_circ')
    % print -f1 -dpsc2 anterior_tension_plot_circ.eps
    print(fig_anterior_circ, '-dpsc2', 'anterior_tension_plot_circ');
    
    figure(fig_anterior_rad); 
    axis(limits_horiz_view);
    set(fig_anterior_rad, 'Position', [100, 100, x_pixels, y_pixels])
    set(fig_anterior_rad, 'PaperPosition', [100, 100, x_pixels, y_pixels])
    %set(fig_anterior_rad,'PaperPositionMode','auto')
    axis off; 
    printfig(fig_anterior_rad, 'anterior_tension_plot_radial')
    % print -f1 -dpsc2 anterior_tension_plot_radial.eps
    % print(fig_anterior_rad, '-dpsc2', 'anterior_tension_plot_radial');
    
    

    % printfig(fig_anterior_both, 'anterior_tension_plot_both')

    figure(fig_posterior_circ); 
    axis(limits_horiz_view); 
    printfig(fig_posterior_circ, 'posterior_tension_plot_circ')

    figure(fig_posterior_rad); 
    axis(limits_horiz_view); 
    printfig(fig_posterior_rad, 'posterior_tension_plot_radial')

    % printfig(fig_posterior_both, 'posterior_tension_plot_both')
    
    
    
%     fig = figure; 
%     
%     n_colors = 100; 
%     colormap(make_colormap(n_colors)); 
% 
%     n_ticks = 5; 
%     tick_array = linspace(0,1,n_ticks); 
%     tick_labels = {}; 
%     for i=1:length(tick_array)
%         tick=tick_array(i); 
%         tension = tick * max_tension * 1e-3; 
%         tick_labels{i} = sprintf('%.2f', tension); 
%     end 
% 
%     cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels); 
%     cbar.Label.String = 'Tension (K Dyne)'; 
%     
%     printfig(fig, 'colorbar_only')

end 




total_tension_plots = true; 
if total_tension_plots

    debug = false; 
    if debug
        n = 32; 
    else
        n = 512;
    end 
    
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
    total_tension_fig_anterior = total_tension_surf_plot(valve.leaflets(1), anterior); 
    axis equal; 
    axis off; 
    printfig(total_tension_fig_anterior, 'total_tension_fig_anterior'); 
    
    anterior = false; 
    total_tension_fig_posterior = total_tension_surf_plot(valve.leaflets(1), anterior); 
    axis equal;
    axis off; 
    printfig(total_tension_fig_posterior, 'total_tension_fig_posterior'); 

end 






