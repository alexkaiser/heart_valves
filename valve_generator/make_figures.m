% 
% script to output many figures for thesis 
% 


N = 32; 

radial       = true; 
bead_slip    = true; 
attached     = false; 
leaflet_only = false; 
optimization = false; 
decreasing_tension = true; 

valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, decreasing_tension); 


output_mesh_schematic(valve); 


