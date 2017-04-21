function [valve] = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential, decreasing_tension)
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
valve.max_it                = 4000; 
valve.max_it_continuation   = 2000; 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 5;  
valve.max_total_fails       = 40; 

if exist('attached', 'var') 
    valve.attached = attached; 
else 
    valve.attached = false; 
end 

% Valve skeleton parameters 

% original, supposedly diastolic but we have been using as systolic 
% valve.r = 1.606587877768772; 
% valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
% valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 

% valve.r = 2.188524100000000; 
% 
% valve.left_papillary  = [-2.307266247008847; -1.497270564906998; -2.639154662183959];
% valve.right_papillary = [-2.058291813251097;  1.497270564906999; -2.712241962241507]; 
% 
% 
% % Places papillary attachments in linear interpolant between single point tips 
% % vector pointing along line from left to right papillary 
% l_to_r_papillary = (valve.right_papillary - valve.left_papillary); 
% l_to_r_papillary = l_to_r_papillary / norm(l_to_r_papillary);
% 
% 
% valve.papillary_radius = 0.25; 
%  
% valve.left_papillary_center  = valve.left_papillary  + valve.papillary_radius * l_to_r_papillary; 
% valve.right_papillary_center = valve.right_papillary - valve.papillary_radius * l_to_r_papillary; 


diastolic_increment = [0; 0; 0.0]; 
valve.left_papillary_diastolic  = []; 
valve.right_papillary_diastolic = []; 


split_papillary = true; 
valve.split_papillary = split_papillary; 
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


valve.decreasing_tension = decreasing_tension; 

if decreasing_tension
 
    % good total value (not including mesh parameters) at N=32
    dec_tension_coeff_32 = 0.002238985441466; 
    
    dec_tension_coeff_base = dec_tension_coeff_32 * 32^2; 
    
    valve.c_dec_tension_circumferential = 2.0 * dec_tension_coeff_base; 
    valve.c_dec_tension_radial          = 2.0 * dec_tension_coeff_base; 
    valve.c_dec_tension_chordae         = 2.0 * dec_tension_coeff_base; 
else 
    valve.dec_tension  = 0.0; 
end 


valve.diff_eqns = @difference_equations_bead_slip; 
valve.jacobian  = @build_jacobian_bead_slip;
        

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

% box width 
valve.L = 2.5; 

MMHG_TO_CGS      = 1333.22368;
valve.p_physical = 110 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

% scaling for target points 
valve.target_multiplier = 40/128; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
valve.n_lagrangian_tracers = 0; 

% Uses configuration of X 
valve.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
valve.num_copies = 1; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_constitutive = false; 

% no reflections in this version 
reflect_x = false; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
% Always true in this version 
radial_and_circumferential = true; 

% physical units create a scalar multiple of the old 
% this multiple is large number, so we want to scale the old tolerance accordingly 
% 8.3326e-04 is a good number here
valve.tol_global = 1e-3;


% controls initial guess tree vertex placement 
tree_frac = 0.5;

% name of structure 
name = 'leaflet'; 

% ring to ring fibers, always zero for now 
ring_to_ring_range = 0; 




sytole_skeleton = false; 
if sytole_skeleton 
    valve.skeleton = valve_points_ct_systole(); 
else 
    valve.skeleton = valve_points_ct_diastole(); 
end 


valve.r        = valve.skeleton.r; 


% Approximate Lagrangian mesh spacing at ring 
% Used for later splitting of springs 
% If any spring is placed at more than double this length an extra vertex is placed
valve.ds = 2*pi*valve.r / N; 


left_papillary_idx  = 1; 
right_papillary_idx = 2; 

valve.dip_anterior_systole = false; 
valve.r_dip = 0.5; 
valve.total_angle_dip = pi; 


% Base constants, individual pieces are tuned relative to these values

% pressure / tension coefficient ratio
% this tension coefficient is the maximum tension that a fiber can support
valve.pressure_tension_ratio = 0.07; % 0.11 * 0.975; 


% base constant for tensions, derived quantity 
valve.tension_base = valve.p_physical / valve.pressure_tension_ratio; 


% Tension coefficients in two directions 
alpha    = 1.0 * valve.tension_base;  % circumferential 
beta     = 1.0 * valve.tension_base;  % radial


% places this many periodic rings above 
n_rings_periodic = 0; %max(1,N/32); 


parameter_values = 4; 


if parameter_values == 1; 

    % explicit commissural leaflets here 

    % Leaflet mesh has irregular bottom edge 
    % 

    % Commissural leaflets centered at -pi/2, pi/2
    N_comm      = 0; % N/8; 

    % Half of each commissural leaflet takes away from the Anterior leaflets half 
    N_anterior  = N/2; 

    % Posterior takes whatever is left 
    N_posterior = N - 2*N_comm - N_anterior; 

    N_per_direction   = (1/2) * [N_anterior, N_anterior, N_comm, N_comm, N_posterior, N_posterior, N_comm, N_comm]; 

    % Anterior goes down then up 
    leaflet_direction = [-1, 1]; 

    % Right commissural goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % Finally, left commissural leaflet goes down to meet initial point 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % offset from N/2 in initial placement 
    leaflet_N_start = 0; %-N_comm/2 + 1; 

    % Anterior leaflet parameters 
    N_anterior = N/2; 

    % Leaf tensions are all modified 
    valve.leaf_tension_base = 0.5 * valve.tension_base; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    valve.root_tension_base = 0.6 * valve.tension_base; 


    n_trees_anterior = 4; 

    k_0_1_anterior = 1.0 * valve.leaf_tension_base / n_trees_anterior; 

    % vector version 
    k_0_1_anterior = k_0_1_anterior * [1.2; 1; 1; 1.2]; 
    k_root_anterior = 0.7 * valve.root_tension_base / n_trees_anterior; 
    k_root_anterior = k_root_anterior * [1.2; 1; 1; 1.2]; 


    papillary_anterior = zeros(3,n_trees_anterior); 

    n_points = n_trees_anterior/2; 

    left_papillary_range = 1:(n_trees_anterior/2); 
    right_papillary_range  = left_papillary_range + (n_trees_anterior/2);

    papillary_anterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points,  0*pi/4,    pi/4); 
    papillary_anterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,   -pi/4, -0*pi/4); 

    n_leaves_anterior  = N_anterior/n_trees_anterior * ones(n_trees_anterior, 1); 
    % leaflet_direction_anterior = [-1; -1; 1; 1];


    N_posterior = N/2; 

    n_trees_posterior = 8; 

    k_0_1_posterior  = 0.2 * valve.leaf_tension_base; 
    k_0_1_posterior  = k_0_1_posterior * ones(n_trees_posterior,1); 
    k_root_posterior = 0.8 * valve.root_tension_base / n_trees_posterior; 
    k_root_posterior = k_root_posterior * [1; 1; 1.2; 1.2; 1.2; 1.2; 1; 1]; 

    papillary_posterior = zeros(3,n_trees_posterior); 

    n_points = n_trees_posterior/2; 

    right_papillary_range = 1:(n_trees_posterior/2); 
    left_papillary_range  = right_papillary_range + (n_trees_posterior/2); 

    papillary_posterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,    pi/4,  5*pi/4); 
    papillary_posterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points, -5*pi/4,   -pi/4);

    % this is generally pretty good 
    n_leaves_posterior = N_posterior/n_trees_posterior * ones(n_trees_posterior, 1); 
    % leaflet_direction_posterior = [-1; 1; -1; -1; 1; 1; -1; 1]; 

    % concatenate all relevant arrays
    n_leaves           = [n_leaves_anterior; n_leaves_posterior];
    papillary          = [papillary_anterior, papillary_posterior]; 
    k_0_1              = [k_0_1_anterior; k_0_1_posterior]; 
    k_root             = [k_root_anterior; k_root_posterior]; 

elseif parameter_values == 1  
        
    % No explicit commissural leaflet here 
    N_anterior = N/2; 

    total_angle_anterior = 5*pi/6; 

    % Posterior takes whatever is left 
    N_posterior = N - N_anterior; 

    N_per_direction   = [N_anterior/2, N_anterior/2, N_posterior/2, N_posterior/2]; 

    % Anterior goes down then up 
    leaflet_direction = [-1, 1]; 

    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % No offset, starting at commissure 
    leaflet_N_start = 0; 


    % Leaf tensions are all modified 
    valve.leaf_tension_base = 0.5 * valve.tension_base; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    valve.root_tension_base = 0.5 * 0.5905 * valve.tension_base; 


    n_trees_anterior = 4; 

    k_0_1_anterior = 0.8 * 2.0 * valve.leaf_tension_base / n_trees_anterior; 

    % vector version 
    k_0_1_anterior  = k_0_1_anterior * [1; 1; 1; 1]; 
    k_root_anterior = 0.9 * 2.0 * valve.root_tension_base / n_trees_anterior; 
    k_root_anterior = k_root_anterior * [1; 1; 1; 1]; 


    papillary_anterior = zeros(3,n_trees_anterior); 

    n_points = n_trees_anterior/2; 

    left_papillary_range = 1:(n_trees_anterior/2); 
    right_papillary_range  = left_papillary_range + (n_trees_anterior/2);

    papillary_anterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points,  0*pi/4,    pi/4); 
    papillary_anterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,   -pi/4, -0*pi/4); 

    n_leaves_anterior  = N_anterior/n_trees_anterior * ones(n_trees_anterior, 1); 

    n_trees_posterior = 8; 

    k_0_1_posterior  = 0.2 * valve.leaf_tension_base; 
    k_0_1_posterior  = k_0_1_posterior * ones(n_trees_posterior,1); 
    k_root_posterior = 0.8 * 2.0 * valve.root_tension_base / n_trees_posterior; 
    k_root_posterior = k_root_posterior * [1; 1; 1; 1; 1; 1; 1; 1]; 

    papillary_posterior = zeros(3,n_trees_posterior); 

    n_points = n_trees_posterior/2; 

    right_papillary_range = 1:(n_trees_posterior/2); 
    left_papillary_range  = right_papillary_range + (n_trees_posterior/2); 

    papillary_posterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,    pi/4,  5*pi/4); 
    papillary_posterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points, -5*pi/4,   -pi/4);

    % this is generally pretty good 
    n_leaves_posterior = N_posterior/n_trees_posterior * ones(n_trees_posterior, 1); 


    % concatenate all relevant arrays
    n_leaves           = [n_leaves_anterior; n_leaves_posterior];
    papillary          = [papillary_anterior, papillary_posterior]; 
    k_0_1              = [k_0_1_anterior; k_0_1_posterior]; 
    k_root             = [k_root_anterior; k_root_posterior]; 



elseif parameter_values == 3 
    
    
    % Commissural leaflets centered at -pi/2, pi/2
    N_comm      = N/4; 

    % Half of each commissural leaflet takes away from the Anterior leaflets half 
    N_anterior  = N/2 - N_comm; 

    % Posterior takes whatever is left 
    N_posterior = N - 2*N_comm - N_anterior; 

    N_per_direction   = [N_comm/4, N_comm/4, ...
                                 N_anterior/2, N_anterior/2, ...
                                 N_comm/4, N_comm/4, N_comm/4, N_comm/4, ...
                                 N_posterior/2, N_posterior/2, ...
                                 N_comm/4, N_comm/4]; 

    % First position is the midpoint of the left commissural leaflet
    % Leaflet moves up from there to the commissure 
    
    % if true, commissural free points change height in mesh 
    % if false commissural points are flat 
    up_and_down_comm = 1;
    
    leaflet_direction = [0 1*up_and_down_comm]; 

    % Anterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % Right commissural goes down then up 
    leaflet_direction = [leaflet_direction, -1*up_and_down_comm, 0, 0, 1*up_and_down_comm]; 

    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % Finally, left commissural leaflet goes down to meet initial point 
    leaflet_direction = [leaflet_direction, -1*up_and_down_comm, 0]; 

    % offset from N/2 in initial placement 
    leaflet_N_start = (-N_comm/4 + 1) * up_and_down_comm; 
    
    
    % tree parameters 
    
    % tip gets a small tree 
    N_tip = N_comm/4; 
    
    % flat commissure starting at center 
    n_leaves = N_comm/4; 
    
    % up commissure and most of anterior are connected 
    n_leaves = [n_leaves, N_anterior/2]; 
    
    % Anterior tip, both sides 
    n_leaves = [n_leaves, N_tip, N_tip]; 
    
    % Anterior and right commissure 
    n_leaves = [n_leaves, N_anterior/2]; 
    
    % flat commissrue 
    n_leaves = [n_leaves, N_comm/4, N_comm/4]; 
    
    % posterior commissure connection, tip, and next commissure 
    n_leaves = [n_leaves, N_posterior/2, N_tip, N_tip, N_posterior/2]; 
    
    % final flat commissure 
    n_leaves = [n_leaves, N_comm/4]; 
    
    
    % four anterior, four posterior, two commissural flat ones  
    total_trees = 12;  
    
    % get all points needed from left and right papillary locations 
    % these are all placed counterclockwise with respect to the entire setup 
    trees_per_side = total_trees/2; 
    papillary_left  = get_papillary_coords(valve, left_papillary_idx,  trees_per_side, -5*pi/4, pi/4); 
    papillary_right = get_papillary_coords(valve, right_papillary_idx, trees_per_side,   -pi/4, 5*pi/4);
    
    % final three left go with the commissure, then the anterior 
    trees_anterior_to_midline_on_left = 3;
    papillary = papillary_left(:, (trees_per_side - trees_anterior_to_midline_on_left + 1): trees_per_side); 
    
    % then all the right trees 
    papillary = [papillary, papillary_right]; 
    
    % then remaining left trees 
    papillary = [papillary, papillary_left(:, 1:(trees_per_side - trees_anterior_to_midline_on_left))]; 
    
    
    
    % Leaf tensions are all modified 
    valve.leaf_tension_base = 0.5 * valve.tension_base; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    valve.root_tension_base = 0.5 * valve.tension_base; 
    
    
    % Anterior leaflet parameters 
    n_trees_anterior = 4; 
    k_0_1_anterior  = 1.0 * valve.leaf_tension_base / n_trees_anterior; 
    k_0_1_anterior  = k_0_1_anterior * [1.1; .9; .9; 1.1]; 
    k_root_anterior = 0.8 * valve.root_tension_base / n_trees_anterior; 
    k_root_anterior = k_root_anterior * [1.2; .8; .8; 1.2]; 

    
    % posterior parameters 
    n_trees_posterior = 4; 
    k_0_1_posterior  = 0.8 * valve.leaf_tension_base / n_trees_posterior; 
    k_0_1_posterior  = k_0_1_posterior * [1.1; .9; .9; 1.1]; 
    k_root_posterior = 0.8 * valve.root_tension_base / n_trees_posterior; 
    k_root_posterior = k_root_posterior * [1; .8; .8; 1]; 
    
    % commissural leaflet parameters 
    n_trees_comm     = 4; 
    k_0_1_comm       = 0.4  * valve.leaf_tension_base / n_trees_comm; 
    k_root_comm      = 0.25 * valve.leaf_tension_base / n_trees_comm; 
    
    k_0_1  = [k_0_1_comm;  k_0_1_anterior;  k_0_1_comm;  k_0_1_comm;  k_0_1_posterior;  k_0_1_comm ]; 
    k_root = [k_root_comm; k_root_anterior; k_root_comm; k_root_comm; k_root_posterior; k_root_comm]; 
    
    
elseif parameter_values == 4 
    
    % No explicit commissural leaflet here 
    N_anterior = N/2; 

    total_angle_anterior = 5*pi/6; 

    % Posterior takes whatever is left 
    N_posterior = N - N_anterior; 

    N_per_direction   = [N_anterior/2, N_anterior/2, N_posterior/2, N_posterior/2]; 

    % Anterior goes down then up 
    leaflet_direction = [-1, 1]; 

    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % No offset, starting at commissure 
    leaflet_N_start = 0; 


    % Leaf tensions are all modified 
    valve.leaf_tension_base = .9 * valve.tension_base; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    valve.root_tension_base = .9 * 0.5905 * valve.tension_base; 


    n_trees_anterior = 2; 

    k_0_1_anterior = valve.leaf_tension_base / n_trees_anterior; 

    % vector version 
    k_0_1_anterior  = k_0_1_anterior * [1; 1]; 
    k_root_anterior = valve.root_tension_base / n_trees_anterior; 
    k_root_anterior = k_root_anterior * [1; 1]; 


    papillary_anterior = zeros(3,n_trees_anterior); 

    n_points = n_trees_anterior/2; 

    left_papillary_range = 1:(n_trees_anterior/2); 
    right_papillary_range  = left_papillary_range + (n_trees_anterior/2);

    papillary_anterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points,  0*pi/4,    pi/4); 
    papillary_anterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,   -pi/4, -0*pi/4); 

    n_leaves_anterior  = N_anterior/n_trees_anterior * ones(n_trees_anterior, 1); 

    n_trees_posterior = 2; 

    k_0_1_posterior  = valve.leaf_tension_base / n_trees_posterior; 
    k_0_1_posterior  = k_0_1_posterior * [1; 1]; 
    k_root_posterior = valve.root_tension_base / n_trees_posterior; 
    k_root_posterior = k_root_posterior * [1; 1]; 

    papillary_posterior = zeros(3,n_trees_posterior); 

    n_points = n_trees_posterior/2; 

    right_papillary_range = 1:(n_trees_posterior/2); 
    left_papillary_range  = right_papillary_range + (n_trees_posterior/2); 

    papillary_posterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points,    pi/4,  5*pi/4); 
    papillary_posterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points, -5*pi/4,   -pi/4);

    % this is generally pretty good 
    n_leaves_posterior = N_posterior/n_trees_posterior * ones(n_trees_posterior, 1); 


    % concatenate all relevant arrays
    n_leaves           = [n_leaves_anterior; n_leaves_posterior];
    papillary          = [papillary_anterior, papillary_posterior]; 
    k_0_1              = [k_0_1_anterior; k_0_1_posterior]; 
    k_root             = [k_root_anterior; k_root_posterior];  
    
    
end 


leaflet = initialize_leaflet_bead_slip(name,                         ... 
                                N,                                   ...
                                reflect_x,                           ... 
                                total_angle_anterior,                ...
                                papillary,                           ... 
                                n_leaves,                            ...
                                leaflet_direction,                   ...
                                leaflet_N_start,                     ...
                                N_per_direction,                     ...
                                radial_and_circumferential,          ...  
                                alpha,                               ... 
                                beta,                                ... 
                                p_0,                                 ... 
                                k_0_1,                               ... 
                                k_root,                              ... 
                                tree_frac,                           ... 
                                leaflet_only,                        ...
                                ring_to_ring_range,                  ...
                                n_rings_periodic,                    ...
                                valve);  

valve.leaflets(1) = leaflet; 
    
% valve_plot(valve); 
% pause(.1); 

disp('Done with initialize.'); 


