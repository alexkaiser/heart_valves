function [valve] = initialize_valve_data_structures_radial_bead_slip(N)
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
valve.radial_and_circumferential = true; 
valve.bead_slip = true; 


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
valve.num_copies = 3; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_springs_leaflet = false; 

extra_posterior = pi/6; 


% posterior leaflet data structure 
reflect_x   = true; 
total_angle_anterior = pi - extra_posterior; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
radial_and_circumferential = true; 


% Spring constants in two directions 
alpha    = 1.0; 
beta     = 1.0; 
p_0      = 0.0; % no pressure for now 
ref_frac = 0.7; % generic spring constants reduced by this much 

% Chordae parameters 
k_0          = 4.0; 
k_multiplier = 2.0; 
tree_frac    = 0.5;

valve.anterior = initialize_leaflet_bead_slip(N,                  ... 
                                    reflect_x,                    ... 
                                    total_angle_anterior,         ...    
                                    valve.r,                      ... 
                                    valve.left_papillary,         ... 
                                    valve.right_papillary,        ... 
                                    radial_and_circumferential,   ...  
                                    alpha,                        ... 
                                    beta,                         ... 
                                    p_0,                          ... 
                                    ref_frac,                     ...  
                                    k_0,                          ... 
                                    k_multiplier,                 ... 
                                    tree_frac);  


                           

valve.posterior = generate_opposite_leaflet(valve.anterior); 

fig = figure; 
valve_plot(valve); 

'here'















