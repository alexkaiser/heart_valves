function [leaflet valve] = initialize_leaflet_bead_slip(name,                         ...
                                                        N,                            ... 
                                                        reflect_x,                    ...  
                                                        angles,                       ... 
                                                        papillary,                    ...
                                                        n_leaves,                     ...
                                                        tree_n_start,                 ... 
                                                        leaflet_direction,            ...
                                                        leaflet_N_start,              ...
                                                        N_per_direction,              ...
                                                        radial_and_circumferential,   ...  
                                                        tension_coeffs,               ... 
                                                        p_0,                          ...
                                                        n_rings_periodic,             ... 
                                                        n_edge_connectors,            ... 
                                                        valve)
%
% Builds leaflet data structures 
% 
% Input:                                   
%                                   
%     N                             Size parameter for leaflet 
%     reflect_x                     Reflects x coordinate in output if true 
%     total_angle                   Angle on valve ring 
%     a                             Cone filter width 
%     h                             Cone filter height 
%     r                             Mitral ring radius  
%     left_papillary                Left papillary point
%     right_papillary               Right papillary point 
%     radial_and_circumferential    Radial and circumferential fibers if true, diagonal if false 
%     alpha                         Spring constant for fibers of constant k, 'u type'
%     beta                          Spring constant for fibers of constant j, 'v type'
%     p_0                           Pressure, if nonzero solve will include this 
%     ref_frac                      Rest lengths uniformly reduced by this much in solve                   
%     k_0                           Base spring constant in leaves of chordae tree 
%     k_multiplier                  Each subsequent 
%     tree_frac                     Determines the weighting for averages in tree initial conditions 
% 
% Output 
% 
%     leaflet                       Fully initialized leaflet data structure 
%                                   
  
leaflet.name               = name; 
leaflet.N                  = N; 

leaflet.num_trees          = size(papillary, 2); 
leaflet.n_rings_periodic   = n_rings_periodic; 
leaflet.n_edge_connectors  = n_edge_connectors; 


leaflet.decreasing_tension = valve.decreasing_tension;
if ~leaflet.decreasing_tension         
    leaflet.c_dec_tension_chordae = 0; 
end 


leaflet.reflect_pressure    = reflect_x; 

leaflet.diff_eqns = @difference_equations_bead_slip; 
leaflet.jacobian  = @build_jacobian_bead_slip;


leaflet.total_angle_anterior = angles.anterior; 

if valve.commissural_leaflets
    leaflet.total_angle_posterior = angles.posterior; 
end 


leaflet.du = 1/N; 

leaflet.r = valve.r; 

if size(papillary, 2) ~= leaflet.num_trees
    error('Must have as many papillary coordinates as trees'); 
end 

if length(n_leaves) ~= leaflet.num_trees
    error('Must have a size for every tree'); 
end 

if length(n_leaves) ~= leaflet.num_trees
    error('Must have a size for every tree'); 
end 


leaflet.n_leaves                  = n_leaves;
leaflet.leaflet_direction         = leaflet_direction; 
leaflet.leaflet_N_start           = leaflet_N_start; 
leaflet.N_per_direction           = N_per_direction; 
leaflet.papillary                 = papillary; 


% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = radial_and_circumferential; 

if ~radial_and_circumferential
    error('slip model only implemented for radial and circumferential')
end 

leaflet = get_free_edge_ranges_bead_slip(leaflet, tree_n_start);

% information about geometry 
leaflet = get_util_arrays_bead_slip(leaflet, valve); 

% build actual data structure 
leaflet.X = build_initial_fibers_bead_slip(leaflet, valve); 

% Scalar pressure to support 
leaflet.p_0 = p_0; 

% Total number of internal leaflet coordinates (three times number of vertices)
leaflet.total_internal_leaflet    = 3*sum(leaflet.is_internal(:)); 

% Running total number of coordinates including trees 
% Updated as trees are added 
leaflet.total_internal_with_trees = 3*sum(leaflet.is_internal(:)); 


% add trees 
leaflet.chordae_tree = true; 
for tree_idx = 1:leaflet.num_trees
    leaflet = add_chordae(leaflet, tree_idx); 
end 
    
% set coefficients on tensions
[leaflet, valve] = set_tension_coeffs(leaflet, valve, tension_coeffs); 

% parameter structure for collagen based nonlinear constitutive 
if valve.collagen_constitutive
    leaflet.collagen_constitutive = true; 
    leaflet.collagen_curve        = get_collagen_curve_parameters(); 
end 


