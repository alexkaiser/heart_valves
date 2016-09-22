function leaflet = initialize_leaflet(N,                            ... 
                                      reflect_x,                    ... 
                                      total_angle,                  ...    
                                      a,                            ... 
                                      h,                            ... 
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
                                      tree_frac,                    ...
                                      trapezoidal_flat_points)
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
leaflet.reflect_x   = reflect_x; 
leaflet.total_angle = total_angle; 
leaflet.min_angle   = -leaflet.total_angle/2.0; 
leaflet.max_angle   =  leaflet.total_angle/2.0; 

leaflet.filter.a = a; 
leaflet.filter.h = h; 
leaflet.filter.r = r; 

if leaflet.reflect_x
    leaflet.left_papillary  = [-1; 1; 1] .* left_papillary; 
    leaflet.right_papillary = [-1; 1; 1] .* right_papillary; 
else 
    leaflet.left_papillary  = left_papillary; 
    leaflet.right_papillary = right_papillary; 
end 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = radial_and_circumferential; 

if exist('trapezoidal_flat_points', 'var')
    if mod(trapezoidal_flat_points,2) == 0
        leaflet.trapezoidal_flat_points    = trapezoidal_flat_points;
    else 
        error('Must add even number of trapezoidal points'); 
    end 
end 

[leaflet.j_max leaflet.k_max leaflet.free_edge_idx_left leaflet.free_edge_idx_right leaflet.chordae_idx_left leaflet.chordae_idx_right] = get_free_edge_ranges(leaflet);

% information about geometry 
[leaflet.is_internal leaflet.is_bc leaflet.linear_idx_offset leaflet.point_idx_with_bc] = get_util_arrays(leaflet); 

% Reference configuration 
leaflet.R = build_reference_surface(leaflet); 

% Initial configuration is reference configuration 
leaflet.X = leaflet.R;  

% Spring constants in two directions 
leaflet.alpha    = alpha; 
leaflet.beta     = beta; 
leaflet.p_0      = p_0; % no pressure for now 
leaflet.ref_frac = ref_frac; % generic spring constants reduced by this much 

% chordae data structures  
leaflet.chordae_tree = true; 
leaflet.k_0          = k_0; 
leaflet.k_multiplier = k_multiplier; 
leaflet.tree_frac    = tree_frac;
leaflet.chordae      = add_chordae(leaflet); 







