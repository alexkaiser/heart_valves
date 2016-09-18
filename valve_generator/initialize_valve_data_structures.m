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
valve.tol_global = 1e-10;
valve.max_it = 40; 


% Valve skeleton parameters 
valve.r = 1.606587877768772; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

% box width 
valve.L = 2.5; 

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
valve.pressure_tension_ratio = 1.5; 


% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
valve.refinement = N/32.0; 

valve.p_physical = 100; 

% scaling for target points 
valve.target_multiplier = 40; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
valve.n_lagrangian_tracers = 8; 

% Uses configuration of X 
valve.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
valve.num_copies = 1; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_springs_leaflet = true; 


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
    [posterior.free_edge_idx_left posterior.free_edge_idx_right posterior.chordae_idx_left posterior.chordae_idx_right] = get_free_edge_ranges(posterior);
else
    error('Radial and circumferential fibers not implemented ')
end 

% information about geometry 
[posterior.is_internal posterior.is_bc posterior.linear_idx_offset] = get_util_arrays(posterior); 

% Reference configuration 
posterior.R = build_reference_surface(posterior); 

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

valve.posterior = posterior; 


% anterior leaflet data structure 
anterior.N           = N; 
anterior.reflect_x   = false; 
anterior.total_angle = pi; 
anterior.min_angle   = -anterior.total_angle/2.0; 
anterior.max_angle   =  anterior.total_angle/2.0; 

anterior.filter.a        = 1.0; 
anterior.filter.h        = 4.0; 
anterior.filter.r        = valve.r; 
anterior.left_papillary  = valve.left_papillary; 
anterior.right_papillary = valve.right_papillary; 


% Radial and circumferential fibers 
% Or diagonally oriented fibers 
anterior.radial_and_circumferential = false; 

if ~anterior.radial_and_circumferential 
    [anterior.free_edge_idx_left anterior.free_edge_idx_right anterior.chordae_idx_left anterior.chordae_idx_right] = get_free_edge_ranges(anterior);
else
    error('Radial and circumferential fibers not implemented ')
end 

% information about geometry 
[anterior.is_internal anterior.is_bc anterior.linear_idx_offset] = get_util_arrays(anterior); 

% Reference configuration 
anterior.R = build_reference_surface(anterior); 

% Initial configuration is reference configuration 
anterior.X = anterior.R;  

% Spring constants in two directions 
anterior.alpha    = 1.0; 
anterior.beta     = 1.0; 
anterior.p_0      = 0.0; % no pressure for now 
anterior.ref_frac = 0.7; % generic spring constants reduced by this much 

anterior.chordae_tree = true; 


if anterior.chordae_tree
    anterior.k_0          = 1.0; 
    anterior.k_multiplier = 1.8;  % 2.0; 
    anterior.tree_frac    = 0.5;
    anterior.chordae      = add_chordae(anterior); 
else 
    error('Non-tree chordae not implemented'); 
end 

valve.anterior = anterior; 














