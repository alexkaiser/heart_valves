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
                                      tree_frac,                    ...
                                      leaflet_only,                 ...
                                      repulsive_potential,          ...
                                      repulsive_power,              ... 
                                      repulsive_coeff)
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
leaflet.total_angle  = total_angle; 
leaflet.leaflet_only = leaflet_only; 


leaflet.repulsive_potential = repulsive_potential; 
leaflet.repulsive_power     = repulsive_power; 
leaflet.repulsive_coeff     = repulsive_coeff; 
leaflet.reflect_pressure    = reflect_x; 

if leaflet_only
    leaflet.diff_eqns = @difference_equations_bead_slip_leaflet_only; 
    leaflet.jacobian  = @build_jacobian_bead_slip_leaflet_only;
else 
    leaflet.diff_eqns = @difference_equations_bead_slip; 
    leaflet.jacobian  = @build_jacobian_bead_slip;
    leaflet.energy    = @energy_bead_slip;
end 

if reflect_x 
    leaflet.min_angle   = pi + leaflet.total_angle/2.0;
    leaflet.max_angle   = pi - leaflet.total_angle/2.0;
else 
    leaflet.min_angle   = -leaflet.total_angle/2.0; 
    leaflet.max_angle   =  leaflet.total_angle/2.0; 
end 

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

[leaflet.j_max leaflet.k_max leaflet.free_edge_idx_left leaflet.free_edge_idx_right leaflet.chordae_idx_left leaflet.chordae_idx_right] = get_free_edge_ranges_bead_slip(leaflet);

% information about geometry 
[leaflet.is_internal leaflet.is_bc leaflet.linear_idx_offset leaflet.point_idx_with_bc] = get_util_arrays_bead_slip(leaflet); 


leaflet.X = build_initial_fibers_bead_slip(leaflet); 

% NaN mask so using bad values will give errors 
leaflet.R = NaN * zeros(size(leaflet.X)); 

free_edge_idx_left  = leaflet.free_edge_idx_left; 
free_edge_idx_right = leaflet.free_edge_idx_right; 

for i=1:size(free_edge_idx_left, 1)
    j = free_edge_idx_left(i,1); 
    k = free_edge_idx_left(i,2); 
    
    % Free edge and neighbors have rest positions 
    % Set to current position for now 
    leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
    leaflet.R(:,j+1,k  ) = leaflet.X(:,j+1,k  );  
    leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);  
end

for i=1:size(free_edge_idx_right, 1)
    j = free_edge_idx_right(i,1); 
    k = free_edge_idx_right(i,2); 
   
    % right free edge has neighbors up in k 
    % but down in j
    leaflet.R(:,j  ,k  ) = leaflet.X(:,j  ,k  );  
    leaflet.R(:,j-1,k  ) = leaflet.X(:,j-1,k  );  
    leaflet.R(:,j  ,k+1) = leaflet.X(:,j  ,k+1);
end 


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

