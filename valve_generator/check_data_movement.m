
N = 4; 

% Main data structure with everything 
valve.N = N; 
valve.tol_global = 1e-10;
valve.max_it = 40; 

% Valve skeleton parameters 
valve.r = 1.606587877768772; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 

% posterior leaflet data structure 
posterior.N           = N; 
posterior.reflect_x   = true; 
posterior.total_angle = pi + pi/6 + pi/12; 
posterior.min_angle   = -posterior.total_angle/2.0; 
posterior.max_angle   =  posterior.total_angle/2.0; 

posterior.filter.a = 1.0; 
posterior.filter.h = 3.0; 
posterior.filter.r = valve.r; 

if posterior.reflect_x
    posterior.left_papillary  = [-1; 1; 1] .* valve.left_papillary; 
    posterior.right_papillary = [-1; 1; 1] .* valve.right_papillary; 
else 
    posterior.left_papillary  = valve.left_papillary; 
    posterior.right_papillary = valve.right_papillary; 
end 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
posterior.radial_and_circumferential = false; 

if ~posterior.radial_and_circumferential 
    [posterior.free_edge_idx_left posterior.free_edge_idx_right] = get_free_edge_ranges(posterior);
else
    error('Radial and circumferential fibers not implemented ')
end 


% Reference configuration 
[posterior.R posterior.is_internal posterior.is_bc] = build_reference_surface(posterior); 

% Initial configuration is reference configuration 
posterior.X = posterior.R;  

% Spring constants in two directions 
posterior.alpha    = 1.0; 
posterior.beta     = 1.0; 
posterior.p_0      = 0.0; % no pressure for now 
posterior.ref_frac = 0.7; % generic spring constants reduced by this much 

posterior.chordae_tree = true; 
if posterior.chordae_tree
    posterior.k_0          = 1.0; 
    posterior.k_multiplier = 1.8;  % 2.0; 
    posterior.tree_frac    = 0.5;
    posterior.chordae      = add_chordae(posterior); 
else 
    error('Non-tree chordae not implemented'); 
end 


'version with no movement, including b.c.s'
posterior.X 


% 'linear order on internal points'
X_linearized = linearize_internal_points(posterior, posterior.X); 


'moved back to 3d ordering'
params = internal_points_to_2d(X_linearized, posterior);  
params.X 


'with chordae before linearization'
posterior.X 
posterior.chordae.C_left
posterior.chordae.C_right 

'linearized with chordae'
X_and_chordae_linearized = linearize_internal_points(posterior, posterior.X, posterior.chordae.C_left, posterior.chordae.C_right) 


'after return to normal data structure'
posterior = internal_points_to_2d(X_and_chordae_linearized, posterior);  
posterior.X 
posterior.chordae.C_left
posterior.chordae.C_right 











