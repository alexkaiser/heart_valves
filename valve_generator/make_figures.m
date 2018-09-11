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
        
    circ = false; 
    radial = true; 
    fig_anterior_rad = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    view(90,0)
    axis equal; 
    anterior_axes_radial = gca
    limits_horiz_view = axis 

    if true % ~debug 
        circ = true; 
        radial = false; 
        fig_anterior_circ = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
        view(90,0)
        axis equal; 
        anterior_axes_circ = gca 
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

        anterior = false; 
        circ = true; 
        radial = false; 
        fig_posterior_circ = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
        view(90,0)
        axis equal;
        posterior_axes_circ = gca 
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

        circ = false; 
        radial = true; 
        fig_posterior_rad = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
        view(90,0)
        axis equal;
        posterior_axes_rad = gca
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 
    end 

    % figure out and manually set the figure output size in pixels 
    lims = axis; 
    x_size = lims(2)-lims(1); 
    y_size = lims(4)-lims(3); 
    x_pixels = floor(100 * x_size); 
    y_pixels = floor(100 * y_size);
    
    figure(fig_anterior_rad); 
    axis(limits_horiz_view);
    %set(fig_anterior_rad, 'Position', [100, 100, x_pixels, y_pixels])
    %set(fig_anterior_rad,'PaperPositionMode','auto')
    axis off; 
    set(fig_anterior_rad, 'Renderer', 'Painters');
    print(fig_anterior_rad, '-depsc', 'anterior_tension_plot_radial_uncropped');
    % export_fig anterior_tension_plot_radial_export_uncropped -eps 
    
    if true % ~debug 

        figure(fig_anterior_circ);
        axis(limits_horiz_view); 
        %set(fig_anterior_circ, 'Position', [100, 100, x_pixels, y_pixels])
        %set(fig_anterior_circ,'PaperPositionMode','auto')
        axis off; 
        set(fig_anterior_circ, 'Renderer', 'Painters');
        print(fig_anterior_circ, '-depsc', 'anterior_tension_plot_circ_uncropped');
        % export_fig anterior_tension_plot_circ_export_uncropped -eps -transparent
    
        figure(fig_posterior_circ); 
        axis(limits_horiz_view);
        %set(fig_posterior_circ, 'Position', [100, 100, x_pixels, y_pixels])
        %set(fig_anterior_rad,'PaperPositionMode','auto')
        axis off 
        set(fig_posterior_circ, 'Renderer', 'Painters');
        printfig(fig_posterior_circ, 'posterior_tension_plot_circ_uncropped')
        % export_fig posterior_tension_plot_circ_export_uncropped -eps -transparent

        figure(fig_posterior_rad); 
        axis(limits_horiz_view);
        %set(fig_posterior_rad, 'Position', [100, 100, x_pixels, y_pixels])
        %set(fig_anterior_rad,'PaperPositionMode','auto')
        axis off
        set(fig_posterior_rad, 'Renderer', 'Painters');
        printfig(fig_posterior_rad, 'posterior_tension_plot_radial_uncropped')
        % export_fig posterior_tension_plot_radial_export_uncropped -eps -transparent 
    end 
    
    
    
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




total_tension_plots = false; 
if total_tension_plots

    debug = false; 
    if debug
        n = 32; 
    else
        n = 512;
    end 
    
    fiber_output    = true; 
    fiber_stride    = 16; 
    stride_offset_j = 4; 
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
    total_tension_fig_anterior = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j); 
    axis equal; 
    axis off; 
    axis tight; 
%     limits_anterior = axis; 
%     lims          = axis; 
    
    pos = get(gcf, 'Position');
    x_pixels_anterior = pos(3); 
    y_pixels_anterior = pos(4); 

    anterior = false; 
    total_tension_fig_posterior = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j); 
    axis equal;
    axis off; 
    axis tight; 
%     limits_posterior = axis; 
%     lims = update_axes(lims, limits_posterior); 
    
    pos = get(gcf, 'Position');
    x_pixels_posterior = pos(3); 
    y_pixels_posterior = pos(4); 
    
%     x_size = lims(2)-lims(1); 
%     y_size = lims(4)-lims(3); 
%     x_pixels = floor(100 * x_size); 
%     y_pixels = floor(100 * y_size);
    
    x_pixels = max(x_pixels_anterior, x_pixels_posterior); 
    y_pixels = max(y_pixels_anterior, y_pixels_posterior); 
    
    fprintf('passed initial calls\n')
    
    figure(total_tension_fig_anterior)
    axis off; 
    set(total_tension_fig_anterior, 'Position', [100, 100, x_pixels, y_pixels])   
    % set(total_tension_fig_anterior,'paperpositionmode','auto')
    % print(total_tension_fig_anterior, '-depsc', '-loose', 'total_tension_fig_anterior'); 
    set(total_tension_fig_anterior, 'Renderer', 'Painters');
    printfig(total_tension_fig_anterior, 'total_tension_fig_anterior_uncropped'); 
    % export_fig total_tension_fig_anterior -eps -transparent
	% export_fig(total_tension_fig_anterior, '-eps', 'total_tension_fig_anterior'); 
    
    fprintf('passed anterior write\n')
    
    figure(total_tension_fig_posterior)
    axis off; 
    set(total_tension_fig_posterior, 'Position', [100, 100, x_pixels, y_pixels])
    set(total_tension_fig_posterior, 'Renderer', 'Painters');
    printfig(total_tension_fig_posterior, 'total_tension_fig_posterior_uncropped'); 
    % export_fig total_tension_fig_posterior -eps -transparent
    fprintf('passed posterior write\n')
    
end 




total_tension_tree_detail = false; 
if total_tension_tree_detail

    debug = false; 
    if debug
        n = 32; 
    else
        n = 512;
    end 
    
    fiber_output    = false; 
    fiber_stride    = 1; 
    stride_offset_j = 0; 
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
    tree_detail_tension_fig_anterior = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j); 
    set(tree_detail_tension_fig_anterior, 'Renderer', 'Painters');
    
    axis equal;
    
    left_view = false;  
    
    if left_view 
        view(-118,20)
        xlabel('x')
        ylabel('y')
        zlabel('z')

        axis([-2.5 0.0 .2 1.4 -2.25 -1.0]);
    else 
        
        % 'pause'
        
        % axis([-2.1 -0.3 -1.7  1.1 -2.8 -1.0])
        axis([-2.5 0.0 -1.4 .2 -2.25 -1.0]);
        
    end
    
    
    
    
    axis off
    
    % saveas(tree_detail_tension_fig_anterior, 'tree_detail_tension_fig_anterior_uncropped.fig'); 
    
%     x_min = -2.5; 
%     x_max = 0; 
%     
%     y_center = .6; 
%     y_radius = .5;
%     
%     y_min = y_center - y_radius; 
%     y_max = y_center + y_radius; 
%     
%     z_center = -2.0; 
%     z_radius = .2; 
%     
%     z_min = z_center - z_radius; 
%     z_max = z_center + z_radius; 
%     
%     axis([x_min x_max y_min y_max z_min z_max]);
%     

%     % zoom(2) 
%     axis off 
    
    output_from_figure = false; 
    if output_from_figure
        tree_detail_tension_fig_anterior = openfig('tree_detail_tension_fig_anterior.fig')
    end
    
    set(tree_detail_tension_fig_anterior, 'Position', [100, 100, 1000, 1000])
    set(tree_detail_tension_fig_anterior,'PaperPositionMode','auto')
    
    if fiber_output
        printfig(tree_detail_tension_fig_anterior, 'tree_detail_tension_fig_anterior_WITH_FIBERS_uncropped'); 
    else 
        printfig(tree_detail_tension_fig_anterior, 'tree_detail_tension_fig_anterior_uncropped'); 
    end 
end 





tension_plots_surf = true; 
if tension_plots_surf

    debug = false; 
    if debug
        n = 32; 
    else
        n = 512;
    end 
    
    fiber_output    = true; 
    fiber_stride    = 16; 
    stride_offset_j = -1; 
    
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
        
    circ = false; 
    radial = true; 
    fig_anterior_rad = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j, circ, radial); 
    view(90,0)
    axis equal; 
    anterior_axes_radial = gca
    limits_horiz_view = axis 

    if ~debug 
        circ = true; 
        radial = false; 
        fig_anterior_circ = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j, circ, radial); 
        view(90,0)
        axis equal; 
        anterior_axes_circ = gca 
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

        anterior = false; 
        circ = true; 
        radial = false; 
        fig_posterior_circ = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j, circ, radial); 
        view(90,0)
        axis equal;
        posterior_axes_circ = gca 
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 

        circ = false; 
        radial = true; 
        fig_posterior_rad = total_tension_surf_plot(valve.leaflets(1), anterior, fiber_output, fiber_stride, stride_offset_j, circ, radial); 
        view(90,0)
        axis equal;
        posterior_axes_rad = gca
        limits_horiz_tmp = axis 
        limits_horiz_view = update_axes(limits_horiz_view, limits_horiz_tmp); 
    end 

    % figure out and manually set the figure output size in pixels 
    lims = axis; 
    x_size = lims(2)-lims(1); 
    y_size = lims(4)-lims(3); 
    x_pixels = floor(100 * x_size); 
    y_pixels = floor(100 * y_size);
    
    figure(fig_anterior_rad); 
    axis(limits_horiz_view);
    axis off; 
    
    set(fig_anterior_rad, 'Renderer', 'Painters');
    print(fig_anterior_rad, '-depsc', 'anterior_tension_plot_radial_surf_uncropped');
    
    if ~debug 

        figure(fig_anterior_circ);
        axis(limits_horiz_view); 
        axis off; 
        set(fig_anterior_circ, 'Renderer', 'Painters');
        print(fig_anterior_circ, '-depsc', 'anterior_tension_plot_circ_surf_uncropped');
    
        figure(fig_posterior_circ); 
        axis(limits_horiz_view);
        axis off 
        set(fig_posterior_circ, 'Renderer', 'Painters');
        printfig(fig_posterior_circ, 'posterior_tension_plot_circ_surf_uncropped')

        figure(fig_posterior_rad); 
        axis(limits_horiz_view);
        axis off
        set(fig_posterior_rad, 'Renderer', 'Painters');
        printfig(fig_posterior_rad, 'posterior_tension_plot_radial_surf_uncropped')
        
    end 

end 

