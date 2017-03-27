function leaflet = initialize_leaflet_bead_slip(N,                  ... 
                                      reflect_x,                    ...  
                                      angles,                       ... 
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
                                      leaflet_only,                 ...
                                      ring_to_ring_range,           ...
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
                
leaflet.N            = N; 
leaflet.leaflet_only = leaflet_only; 

leaflet.tension_base = valve.tension_base; 

leaflet.ring_to_ring_range = ring_to_ring_range; 


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
leaflet.dv = 1/N; 


leaflet.r = r; 

leaflet.left_papillary  = left_papillary; 
leaflet.right_papillary = right_papillary; 


% Radial and circumferential fibers 
% Or diagonally oriented fibers 
leaflet.radial_and_circumferential = radial_and_circumferential; 

if ~radial_and_circumferential
    error('slip model only implemented for radial and circumferential')
end 

leaflet = get_free_edge_ranges_bead_slip(leaflet);

% information about geometry 
leaflet = get_util_arrays_bead_slip(leaflet); 


leaflet.X = build_initial_fibers_bead_slip(leaflet); 

% % NaN mask so using bad values will give errors 
% leaflet.R = NaN * zeros(size(leaflet.X)); 
% 
% free_edge_idx_left  = leaflet.free_edge_idx_left; 
% free_edge_idx_right = leaflet.free_edge_idx_right; 
% 
% for i=1:size(free_edge_idx_left, 1)
%     j = free_edge_idx_left(i,1); 
%     k = free_edge_idx_left(i,2); 
%     
%     % Free edge and neighbors have rest positions 
%     % Set to current position for now 
%     leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
%     leaflet.R(:,j+1,k  ) = leaflet.X(:,j+1,k  );  
%     leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);  
% end
% 
% for i=1:size(free_edge_idx_right, 1)
%     j = free_edge_idx_right(i,1); 
%     k = free_edge_idx_right(i,2); 
%    
%     % right free edge has neighbors up in k 
%     % but down in j
%     leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
%     leaflet.R(:,j-1,k  ) = leaflet.X(:,j-1,k  );  
%     leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);
% end 


% Spring constants in two directions 
leaflet.alpha    = alpha; 
leaflet.beta     = beta; 
leaflet.p_0      = p_0; 
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

