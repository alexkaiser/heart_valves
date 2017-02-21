function [valve] = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential)
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
valve.tol_global = 1e-9;
valve.max_it = 4000; 

if exist('attached', 'var') 
    valve.attached = attached; 
else 
    valve.attached = false; 
end 

% Valve skeleton parameters 
valve.r = 1.606587877768772; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 
valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.leaflet_only = leaflet_only; 
valve.optimization = optimization; 
valve.repulsive_potential = repulsive_potential; 
valve.repulsive_power     = 1; 

if repulsive_potential
 
    % good total value (not including mesh parameters) at N=32
    repulsive_coeff_32 = 0.002238985441466; 
    
    repulsive_coeff_base = repulsive_coeff_32 * 32^2; 
    
    % scale so that when multiplied by above value gives the correct value 
    % valve.repulsive_coeff = repulsive_coeff_32 * 32^2; 
    
    valve.c_repulsive_circumferential = 2.0 * repulsive_coeff_base; 
    valve.c_repulsive_radial          = 6.0 * repulsive_coeff_base; 
    valve.c_repulsive_chordae         = 1.0 * repulsive_coeff_base; 
else 
    valve.repulsive_coeff  = 0.0; 
end 


% most interesting power 2 at this point 
% valve.repulsive_power     = 2; 
% valve.repulsive_coeff     = (0.28) * 1.0e-3;

% function pointers 
if attached 
    valve.diff_eqns = @difference_equations_bead_slip_attached; 
    valve.jacobian  = @build_jacobian_bead_slip_attached; 
    
    if leaflet_only
        error('leaflet_only not implemented for attached')
    end 
        
else 
    if leaflet_only
        valve.diff_eqns = @difference_equations_bead_slip_leaflet_only; 
        valve.jacobian  = @build_jacobian_bead_slip_leaflet_only;
    else 
        valve.diff_eqns = @difference_equations_bead_slip; 
        valve.jacobian  = @build_jacobian_bead_slip;
        
    end 
end 

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

% box width 
valve.L = 2.5; 

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
valve.pressure_tension_ratio = .75; 


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


% anterior leaflet data structure 
reflect_x = false; 
total_angle_anterior = pi - extra_posterior; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
radial_and_circumferential = true; 


% Spring constants in two directions 
alpha    =  1.0;  % circumferential 
beta     =  1.0;  % radial 
p_0      =  0.0;  %-0.12; % negative sign on anterior leaflet 
ref_frac =  0.7;  % generic spring constants reduced by this much 


% Add energy function for zero pressure case 
if (p_0 == 0.0) && (~leaflet_only)
    valve.energy = @energy_bead_slip; 
end 


% Chordae parameters 

% tree has half as many leaves as total number of radial fibers N
N_tree = N/2; 

% base constant for force scaling
% this is the total force, in current units, 
% in the leaf generation of the chordae tree 
% this is an arbitrary constant determined by guess and check 
k_0_1 = 1.8 * 0.5; 

% force on each leaf in the chordae tree 
k_0   = k_0_1 / N_tree; 

% constant tension at the root of the tree 
% this is determined by hand tuning k_multiplier at coarse resolution 
% then taking the k_root 

% base good value 
% k_root = 1.889568000000001e+01 / 32; 

% adjust accordingly
k_root = 0.9 * 1.889568000000001e+01 / 32; 

% multiplier necessary to maintain constant root tension 
% and constant total leaf tension 
k_multiplier = 2.0 * (k_root/k_0_1)^(1/log2(N_tree)); 


tree_frac    = 0.5;

% double tree strength in attached version 
if attached 
    k_0 = k_0 * 2.0; 
end 

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
                                    tree_frac,                    ... 
                                    leaflet_only,                 ...
                                    valve);  

                   

if valve.attached 
    valve.posterior = generate_opposite_leaflet(valve.anterior); 
    fig = figure; 
    
else 
    
    total_angle_posterior = 2*pi - total_angle_anterior; 
    reflect_x = true; 
    
    % reflect pressure also 
    if reflect_x
        p_0 = -p_0; 
    end 
    
    valve.posterior = initialize_leaflet_bead_slip(N,             ...
                                    reflect_x,                    ... 
                                    total_angle_posterior,        ...    
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
                                    tree_frac,                    ... 
                                    leaflet_only,                 ...
                                    valve);  

    
    
end 

% valve_plot(valve); 

'done with initialize'















