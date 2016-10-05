function leaflet = initialize_leaflet_bead_slip(N,                  ... 
                                      reflect_x,                    ...  
                                      total_angle,                  ... 
                                      r,                            ... 
                                      left_papillary,               ... 
                                      right_papillary,              ... 
                                      radial_and_circumferential,   ...  
                                      alpha,                        ... 
                                      beta,                         ... 
                                      p_0,                          ...
                                      ref_frac,                     ...  
                                      k_0,                          ... 
                                      k_multiplier,                 ... 
                                      tree_frac)
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
                                  
leaflet.N           = N; 
leaflet.total_angle = total_angle; 
leaflet.min_angle   = -leaflet.total_angle/2.0; 
leaflet.max_angle   =  leaflet.total_angle/2.0; 

leaflet.r = r; 

leaflet.left_papillary  = left_papillary; 
leaflet.right_papillary = right_papillary; 


% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = radial_and_circumferential; 

if ~radial_and_circumferential
    error('slip model only implemented for radial and circumferential')
end 

[leaflet.j_max leaflet.k_max leaflet.free_edge_idx_left leaflet.free_edge_idx_right leaflet.chordae_idx_left leaflet.chordae_idx_right] = get_free_edge_ranges_bead_slip(leaflet);

% information about geometry 
[leaflet.is_internal leaflet.is_bc leaflet.linear_idx_offset leaflet.point_idx_with_bc] = get_util_arrays_bead_slip(leaflet); 


leaflet.X = build_initial_fibers_bead_slip(leaflet); 


% Spring constants in two directions 
leaflet.alpha    = alpha; 
leaflet.beta     = beta; 
leaflet.p_0      = p_0; % no pressure for now 
leaflet.ref_frac = ref_frac; % generic spring constants reduced by this much 

% chordae data structures  
if exist('k_0', 'var') && exist('k_multiplier', 'var') && exist('tree_frac', 'var')
    leaflet.chordae_tree = true; 
    leaflet.k_0          = k_0; 
    leaflet.k_multiplier = k_multiplier; 
    leaflet.tree_frac    = tree_frac;
    leaflet.chordae      = add_chordae(leaflet); 
else 
    leaflet.chordae_tree = false; 
end 

