function [fig_circ fig_rad fig_ratio] = make_aortic_plots(leaflet, fig_circ, fig_rad, fig_ratio)

if ~exist('fig_circ', 'var')
    fig_circ = figure; 
end 
if ~exist('fig_rad', 'var')
    fig_rad = figure; 
end 
if ~exist('fig_ratio', 'var')
    fig_ratio = figure; 
end 

fiber_output    = true; 
fiber_stride    = 4; 
stride_offset_j = 0; 

height_plot = false; 

circ  = true; 
rad   = false; 
ratio = false; 

set(0, 'CurrentFigure', fig_circ)
[az el] = view; 
clf(fig_circ); 
total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig_circ)
view(az,el);
title('circ tension')

circ  = false; 
rad   = true; 
ratio = false; 
set(0, 'CurrentFigure', fig_rad)
[az el] = view; 
clf(fig_rad); 
total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig_rad)
view(az,el);
title('radial tension')
    
circ  = false; 
rad   = false; 
ratio = true; 
set(0, 'CurrentFigure', fig_ratio)
[az el] = view; 
clf(fig_ratio); 
total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig_ratio)
view(az,el);
title('ratio circ/radial tension')

