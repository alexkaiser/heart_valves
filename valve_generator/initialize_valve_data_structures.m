function [valve] = initialize_valve_data_structures(N)
% 
% Initializes data structures for full solve.  
% 
% Parameters are declared here.
% Should be a script, but want to return in the structures 
% 
% Input: 
%     N   Size parameter used throughout 
% 


% Main data structure with everything 
valve.N = N; 

% Valve skeleton parameters 
valve.r = 1.606587877768772; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 


% posterior leaflet data structure 
valve.posterior.N = N; 
valve.posterior.reflect_z = true; 
valve.posterior.total_angle = pi + pi/6 + pi/12; 
valve.posterior.min_angle   = -valve.posterior.total_angle/2.0; 
valve.posterior.max_angle   =  valve.posterior.total_angle/2.0; 

valve.posterior.filter.a = 1.0; 
valve.posterior.filter.h = 3.0; 
valve.posterior.filter.r = valve.r; 
% valve.posterior.filter.N = N; 
% valve.posterior.filter.min_angle = valve.posterior.min_angle; 
% valve.posterior.filter.max_angle = valve.posterior.max_angle; 

if valve.posterior.reflect_z
    valve.posterior.left_papillary  = [-1; 1; 1] .* valve.left_papillary; 
    valve.posterior.right_papillary = [-1; 1; 1] .* valve.right_papillary; 
else 
    valve.posterior.left_papillary  = valve.left_papillary; 
    valve.posterior.right_papillary = valve.right_papillary; 
end 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
valve.posterior.radial_and_circumferential = false; 

if ~valve.posterior.radial_and_circumferential 
    [valve.posterior.free_edge_idx_left valve.posterior.free_edge_idx_right] = get_free_edge_ranges(valve.posterior);
else
    error('Radial and circumferential fibers not implemented ')
end 


% Reference configuration 
[valve.posterior.R valve.posterior.is_internal valve.posterior.is_bc] = build_reference_surface(valve.posterior); 

% Initial configuration is reference configuration 
valve.posterior.X = valve.posterior.R;  

% Spring constants in two directions 
valve.posterior.alpha    = 1.0; 
valve.posterior.beta     = 1.0; 
valve.posterior.p_0      = 0.0; % no pressure for now 
valve.posterior.ref_frac = 0.7; % generic spring constants reduced by this much 

valve.posterior.chordae_tree = true; 


if valve.posterior.chordae_tree
    valve.posterior.k_0          = 1.0; 
    valve.posterior.k_multiplier = 1.8;  % 2.0; 
    valve.posterior.tree_frac    = 0.5;
    valve.posterior.chordae      = add_chordae(valve.posterior); 
else 
    error('Non-tree chordae not implemented'); 
end 

fig = surf_plot(valve.posterior); 
title('Reference configuration of posterior surface'); 


% anterior leaflet data structure 
% valve.anterior  = struct; 

















