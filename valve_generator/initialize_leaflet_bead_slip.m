function leaflet = initialize_leaflet_bead_slip(name,               ...
                                      N,                            ... 
                                      reflect_x,                    ...  
                                      angles,                       ... 
                                      r,                            ... 
                                      papillary,                    ...
                                      n_leaves,                     ...
                                      tree_direction,               ...
                                      left_papillary_diastolic,     ...
                                      right_papillary_diastolic,    ...
                                      radial_and_circumferential,   ...  
                                      alpha,                        ... 
                                      beta,                         ... 
                                      p_0,                          ...  
                                      k_0_1,                        ... 
                                      k_root,                       ... 
                                      tree_frac,                    ...
                                      leaflet_only,                 ...
                                      ring_to_ring_range,           ...
                                      n_rings_periodic,             ... 
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
leaflet.leaflet_only       = leaflet_only; 
leaflet.tension_base       = valve.tension_base; 
leaflet.ring_to_ring_range = ring_to_ring_range; 

leaflet.num_trees          = size(papillary, 2); 
leaflet.n_rings_periodic   = n_rings_periodic; 


leaflet.repulsive_potential         = valve.repulsive_potential; 
if leaflet.repulsive_potential         
    leaflet.repulsive_power             = valve.repulsive_power; 
    leaflet.c_repulsive_circumferential = valve.c_repulsive_circumferential; 
    leaflet.c_repulsive_radial          = valve.c_repulsive_radial; 
    leaflet.c_repulsive_chordae         = valve.c_repulsive_chordae; 
else 
    leaflet.repulsive_power             = 1; 
    leaflet.c_repulsive_circumferential = 0; 
    leaflet.c_repulsive_radial          = 0; 
    leaflet.c_repulsive_chordae         = 0; 
end


leaflet.decreasing_tension         = valve.decreasing_tension;
if leaflet.decreasing_tension         
    leaflet.c_dec_tension_circumferential = valve.c_dec_tension_circumferential; 
    leaflet.c_dec_tension_radial          = valve.c_dec_tension_radial; 
    leaflet.c_dec_tension_chordae         = valve.c_dec_tension_chordae; 
else 
    leaflet.c_dec_tension_circumferential = 0; 
    leaflet.c_dec_tension_radial          = 0; 
    leaflet.c_dec_tension_chordae         = 0; 
end 


leaflet.reflect_pressure    = reflect_x; 

if leaflet_only
    leaflet.diff_eqns = @difference_equations_bead_slip_leaflet_only; 
    leaflet.jacobian  = @build_jacobian_bead_slip_leaflet_only;
else 
    leaflet.diff_eqns = @difference_equations_bead_slip; 
    leaflet.jacobian  = @build_jacobian_bead_slip;
    leaflet.energy    = @energy_bead_slip;
end 


if length(angles) ~= 2
    error('Must provide min and max angle'); 
end 

leaflet.min_angle = angles(1); 
leaflet.max_angle = angles(2); 

leaflet.du = 1/N; 

leaflet.r = r; 

if size(papillary, 2) ~= leaflet.num_trees
    error('Must have as many papillary coordinates as trees'); 
end 

if length(n_leaves) ~= leaflet.num_trees
    error('Must have a size for every tree'); 
end 

if length(tree_direction) ~= leaflet.num_trees
    error('Must have a size for every tree'); 
end 

if sum(n_leaves) ~= N
    error('Must have a free edge connection for all points on the tree'); 
end 

if sum(tree_direction .* n_leaves) ~= 0
    error('Total leaflet must return to same vertical index from which it started'); 
end 
    


leaflet.n_leaves                  = n_leaves;
leaflet.tree_direction            = tree_direction; 
leaflet.papillary                 = papillary; 

% FIX ME
leaflet.left_papillary_diastolic  = left_papillary_diastolic; 
leaflet.right_papillary_diastolic = right_papillary_diastolic;


% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = radial_and_circumferential; 

if ~radial_and_circumferential
    error('slip model only implemented for radial and circumferential')
end 

leaflet = get_free_edge_ranges_bead_slip(leaflet);

% information about geometry 
leaflet = get_util_arrays_bead_slip(leaflet); 

% layout on 
leaflet.X = build_initial_fibers_bead_slip(leaflet); 


% Spring constants in two directions 
leaflet.alpha    = alpha; 
leaflet.beta     = beta; 
leaflet.p_0      = p_0; 

% Total number of internal leaflet coordinates (three times number of vertices)
leaflet.total_internal_leaflet    = 3*sum(leaflet.is_internal(:)); 

% Running total number of coordinates including trees 
% Updated as trees are added 
leaflet.total_internal_with_trees = 3*sum(leaflet.is_internal(:)); 



leaflet.chordae_tree = true; 
leaflet.k_0_1        = k_0_1; 
leaflet.k_root       = k_root; 
leaflet.tree_frac    = tree_frac;

for tree_idx = 1:leaflet.num_trees
    leaflet = add_chordae(leaflet, tree_idx); 
end 
    


% parameter structure for collagen based nonlinear constitutive 
if valve.collagen_constitutive
    leaflet.collagen_constitutive = true; 
    leaflet.collagen_curve        = get_collagen_curve_parameters(); 
end 


