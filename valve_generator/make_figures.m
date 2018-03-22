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



tension_plots = true; 
if tension_plots

    n = 32;
    
    base_dir  = '/Users/alex/mitral_fully_discrete/valve_generator/meshes/plot_meshes/two_leaflet_8_connector_b7a6aed/'
    file_name = ['mitral_tree_', int2str(n), '_final_data.mat']
    
    load([base_dir, file_name])
  
    anterior = true; 
    circ = true; 
    radial = false; 
    fig_anterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_anterior, 'anterior_tension_plot_circ')
    
    input('Press any key to move on');  
    
    % close(fig_anterior)


    
    circ = false; 
    radial = true; 
    fig_anterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_anterior, 'anterior_tension_plot_radial')
    % close(fig_anterior)

    circ = true; 
    radial = true; 
    fig_anterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_anterior, 'anterior_tension_plot_both')
    % close(fig_anterior)

    anterior = false; 
    circ = true; 
    radial = false; 
    fig_posterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_posterior, 'posterior_tension_plot_circ')
    % close(fig_posterior)

    circ = false; 
    radial = true; 
    fig_posterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_posterior, 'posterior_tension_plot_radial')
    % close(fig_posterior)

    circ = true; 
    radial = true; 
    fig_posterior = fiber_tension_surf_plot(valve.leaflets(1), anterior, circ, radial); 
    savefig(fig_posterior, 'posterior_tension_plot_both')
    % close(fig_posterior)

end 



