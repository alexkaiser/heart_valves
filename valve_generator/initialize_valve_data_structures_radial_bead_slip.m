function [valve] = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, decreasing_tension)
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

% effective infinity by default 
valve.max_it                = 1e8; 
valve.max_continuations     = 1e8; 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 0;  
valve.max_total_fails       = 0; 

if exist('attached', 'var') 
    valve.attached = attached; 
else 
    valve.attached = false; 
end 

split_papillary = true; 
valve.split_papillary = split_papillary; 
valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.leaflet_only = leaflet_only; 
valve.optimization = optimization; 

valve.decreasing_tension = decreasing_tension; 


valve.diff_eqns = @difference_equations_bead_slip; 
valve.jacobian  = @build_jacobian_bead_slip;
        

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

MMHG_TO_CGS      = 1333.22368;
valve.p_physical = 100 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 


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

% add flags to spring files 
% to view and output with a stride 


valve.output.leaflets       = [1;0;1]; 
valve.output.stride_leaflet = N/32; 
valve.output.chordae        = [1;1;1]; 
valve.output.mesh           = [1;0;0]; 
valve.output.cartesian_mesh = [0;0;0]; 
valve.output.stride_mesh    = N/32; 

valve.papillary_movement_times = [0 .1 .46 .5 .8]; 


% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_constitutive = true; 

% Constant strain of pressurized configuration 
valve.strain = .16; 

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

% name of structure 
name = 'leaflet'; 

left_papillary_idx  = 1; 
right_papillary_idx = 2; 


explicit_comm_leaflets = false; 
    

if ~explicit_comm_leaflets 

    % commissural tree version 
    % but without explicit commissural leaflets 
    
    valve.dip_anterior_systole = true; 
    valve.r_dip = 0.75; 
    valve.total_angle_dip = pi; 

    valve.L = 3.0; 
    
    low_papillary = false; 
    tip_radius = .2; 
    valve.skeleton = valve_points_ct_systole(low_papillary, tip_radius); 
    
    
    valve.diastolic_increment = [1.25; 0.0; 0.25]; 

    
    zero_radius = false; 
    if zero_radius
        for i = 1:length(valve.skeleton.papillary)
            valve.skeleton.papillary(i).radius = 0; 
        end 
    end 
    
    vertical_normal_papillary = true; 
    if vertical_normal_papillary 
        for i = 1:length(valve.skeleton.papillary)
            valve.skeleton.papillary(i).normal = [0; 0; 1]; 
        end 
    end 

    % Base constants, individual pieces are tuned relative to these values

    if decreasing_tension
        dec_tension_coeff_base       = 4.6;  
        valve.dec_tension_coeff_base = dec_tension_coeff_base; 
    else 
        valve.dec_tension  = 0.0; 
    end 
    

    
    % tension coefficients structure 

    % pressure / tension coefficient ratio
    % this tension coefficient is the maximum tension that a fiber can support
    % valve.pressure_tension_ratio = 0.055; % 0.11 * 0.975; 
    tension_coeffs.pressure_tension_ratio = 0.055; 
    
    tension_coeffs.dec_tension_coeff_base = 4.6; 

    
    % max tensions in leaflets 
    tension_coeffs.alpha_anterior       = 1.75;  % circumferential 
    tension_coeffs.beta_anterior        = 1.1;  % radial
    tension_coeffs.alpha_posterior      = 1.0;  % circumferential 
    tension_coeffs.beta_posterior       = 1.0;  % radial
    tension_coeffs.alpha_hoops          = 0.5;  % circumferential hoops     
    tension_coeffs.alpha_edge_connector = 1.0;  % circumferential free edge connector 
    tension_coeffs.beta_edge_connector  = 0.01;  % circumferential free edge connector


    % decreasing tension coefficients 
    tension_coeffs.c_circ_dec_anterior       = 2.5;  % circumferential 
    tension_coeffs.c_rad_dec_anterior        = 1.5;  % radial
    tension_coeffs.c_circ_dec_posterior      = 1.0;  % circumferential 
    tension_coeffs.c_rad_dec_posterior       = 1.5;  % radial
    tension_coeffs.c_circ_dec_hoops          = 2.0;  % circumferential hoops
    tension_coeffs.c_rad_dec_hoops_anterior  = 0.5;  % radial hoops, anterior part 
    tension_coeffs.c_rad_dec_hoops_posterior = 0.5;  % radial hoops, posterior part 
    tension_coeffs.c_dec_tension_chordae     = 1.0;  % chordae

    tension_coeffs.c_circ_dec_edge_connector  = 2.0;  % circumferential hoops
    tension_coeffs.c_rad_dec_edge_connector   = 2.0;  % circumferential hoops

    % places this many periodic rings above 
    n_rings_periodic = max(2,N/64); 

    % places circumferential fibers this many below hoops 
    % if the location is not already covered by leaflet
    n_edge_connectors = max(1,N/64);  

    % No explicit commissural leaflet here 
    N_anterior = N/2; 

    angles.anterior = 5*pi/6; 

    % Posterior takes whatever is left 
    N_posterior = N - N_anterior; 

    % store these 
    valve.N_anterior   = N_anterior; 
    valve.N_posterior  = N_posterior;
    valve.commissural_leaflets = false; 
    valve.N_commissure = 0; 


    N_per_direction   = [N_anterior/2, N_anterior/2, N_posterior/2, N_posterior/2]; 

    % Anterior goes down then up 
    leaflet_direction = [-1, 1]; 

    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 

    % No offset, starting at commissure 
    leaflet_N_start = 0; 

    % changes entire tree strength by constants 
    tension_coeffs.tree_tension_multiplier = 1.0; 

    % Leaf tensions are all modified 
    tension_coeffs.leaf_tension_base = .9; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    tension_coeffs.root_tension_base = .9 * 0.5905; 


    n_trees_anterior = 2; 

%     k_0_1_anterior  = 1.1 / n_trees_anterior; 
%     k_0_1_anterior  = k_0_1_anterior * [1; 1]; 
%     k_root_anterior = 1.1 / n_trees_anterior; 
%     k_root_anterior = k_root_anterior * [1; 1]; 

%    n_leaves_anterior  = N_anterior/n_trees_anterior * ones(n_trees_anterior, 1); 

    
    % posterior and included commissural trees 
    n_trees_posterior_and_comm  = 6;
%    n_trees_posterior           = 2; 
%    n_trees_commissure          = 4; 
%    n_trees_commissure_per_side = n_trees_commissure/2; 
%     
%     % include commissural trees in posterior leaflet 
%     n_posterior_tree_total   = N_posterior / 2; 
%     n_commissural_tree_total = N_posterior / 2;
%     
%     n_tree_posterior         = n_posterior_tree_total   / n_trees_posterior; 
%     n_tree_commissure        = n_commissural_tree_total / n_trees_commissure; 
%     
%     k_0_1_posterior          = 0.4 / n_trees_posterior; 
%     k_root_posterior         = 0.4 / n_trees_posterior; 
%     
%     k_0_1_commissure         = 0.5 / n_trees_commissure; 
%     k_root_commissure        = 0.5 / n_trees_commissure; 
% 
%     
%     k_0_1_posterior_and_comm    = zeros(n_trees_posterior_and_comm, 1);
%     k_root_posterior_and_comm   = zeros(n_trees_posterior_and_comm, 1);
%     n_leaves_posterior_and_comm = zeros(n_trees_posterior_and_comm, 1); 
%     
%     
%     j = 1; 
%     for tmp=1:n_trees_commissure_per_side
%         k_0_1_posterior_and_comm(j)    = k_0_1_commissure; 
%         k_root_posterior_and_comm(j)   = k_root_commissure; 
%         n_leaves_posterior_and_comm(j) = n_tree_commissure; 
%         j = j+1; 
%     end 
%     
%     for tmp=1:n_trees_posterior
%         k_0_1_posterior_and_comm(j)    = k_0_1_posterior; 
%         k_root_posterior_and_comm(j)   = k_root_posterior; 
%         n_leaves_posterior_and_comm(j) = n_tree_posterior; 
%         j = j+1; 
%     end    
%     
%     for tmp=1:n_trees_commissure_per_side
%         k_0_1_posterior_and_comm(j)    = k_0_1_commissure; 
%         k_root_posterior_and_comm(j)   = k_root_commissure; 
%         n_leaves_posterior_and_comm(j) = n_tree_commissure; 
%         j = j+1; 
%     end
    
    % this array determines the fraction of N_orig which each tree takes up 
    % this allows us to determine initial fractions of constants that go to each tree 
    frac_of_n_orig = [1/ 4; 1/ 4;  ...   % anterior  
                      1/16; 1/16;  ...   % comm
                      1/ 8; 1/ 8;  ...   % posterior
                      1/16; 1/16];       % comm
    
    % change these to manipulate individial tree coefficients 
    % for sanity reasons, these shuold mostly be one unless you have a good reason to change 
    % note that these are scaled by the fraction of the leaflet that they take up 
    k_0_1_coeff    = frac_of_n_orig .*    ... 
                     [2.2; 2.2;           ...       % anterior  
                      2.0; 2.0;           ...       % anterior and comm, comm and posterior       
                      1.6; 1.6;           ...       % posterior
                      2.0; 2.0];                    % posterior and comm, comm and anterior
                  
    k_root_coeff   = frac_of_n_orig .*    ... 
                    [ 2.2; 2.2;           ...       % anterior  
                      2.0; 2.0;           ...       % anterior and comm, comm and posterior       
                      1.6; 1.6;           ...       % posterior
                      2.0; 2.0];                    % posterior and comm, comm and anterior
                  
                                         
    tension_coeffs.c_dec_chordae_leaf = (1/N)  * [1.0; 1.0; ...       % anterior  
                                                  1.0; 1.0; ...       % comm
                                                  1.0; 1.0; ...       % posterior
                                                  1.0; 1.0];          % comm 
                                              
    % root constants do not scale, because the root 
    % should maintain a consistent length when mesh is changed 
    tension_coeffs.c_dec_chordae_root = (1/256) * [1.0; 1.0; ...       % anterior  
                                                   1.0; 1.0; ...       % comm
                                                   1.0; 1.0; ...       % posterior
                                                   1.0; 1.0];          % comm 
                                              
    
    n_leaves = N * frac_of_n_orig; 
    
    papillary_anterior = zeros(3,n_trees_anterior); 
    n_points = n_trees_anterior/2; 
    left_papillary_range = 1:(n_trees_anterior/2); 
    right_papillary_range  = left_papillary_range + (n_trees_anterior/2);
    
    % arrangements of connection to papillary muscle 
    % angles are measured form approximate x direction on left papillary 
    % negatives and swapped on right papillary 
    min_papillary_angle_anterior = 0; 
    max_papillary_angle_anterior = 0; 
    
    papillary_anterior(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points,  min_papillary_angle_anterior,  max_papillary_angle_anterior); 
    papillary_anterior(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points, -max_papillary_angle_anterior, -min_papillary_angle_anterior); 

    
    papillary_posterior_and_comm = zeros(3,n_trees_posterior_and_comm); 
    n_points = n_trees_posterior_and_comm/2; 
    right_papillary_range = 1:(n_trees_posterior_and_comm/2); 
    left_papillary_range  = right_papillary_range + (n_trees_posterior_and_comm/2); 
    
    % arrangements of anchor points 
    min_papillary_angle_posterior = -pi; 
    max_papillary_angle_posterior = 0; 
    
    papillary_posterior_and_comm(:,right_papillary_range) = get_papillary_coords(valve, right_papillary_idx, n_points, -max_papillary_angle_posterior, -min_papillary_angle_posterior); 
    papillary_posterior_and_comm(:,left_papillary_range)  = get_papillary_coords(valve, left_papillary_idx,  n_points,  min_papillary_angle_posterior,  max_papillary_angle_posterior);

    
    % concatenate all relevant arrays
    % n_leaves           = [n_leaves_anterior; n_leaves_posterior_and_comm];
    papillary          = [papillary_anterior, papillary_posterior_and_comm]; 
    % k_0_1              = [k_0_1_anterior; k_0_1_posterior_and_comm]; 
    % k_root             = [k_root_anterior; k_root_posterior_and_comm]; 
    
    tension_coeffs.k_0_1  = k_0_1_coeff; 
    tension_coeffs.k_root = k_root_coeff; 
    
elseif explicit_comm_leaflets 

    % commissural tree version 
    % with explicit commissural leaflets 
    
    valve.dip_anterior_systole = true; 
    valve.r_dip = 0.75; 
    valve.total_angle_dip = pi; 

    valve.L = 3.0; 
    
    low_papillary = true; 
    tip_radius = .2; 
    valve.skeleton = valve_points_ct_systole(low_papillary, tip_radius); 
    
    
    valve.diastolic_increment = [1.5; 0.0; 0.5]; 

    
    zero_radius = false; 
    if zero_radius
        for i = 1:length(valve.skeleton.papillary)
            valve.skeleton.papillary(i).radius = 0; 
        end 
    end 
    
    vertical_normal_papillary = true; 
    if vertical_normal_papillary 
        for i = 1:length(valve.skeleton.papillary)
            valve.skeleton.papillary(i).normal = [0; 0; 1]; 
        end 
    end 

    % Base constants, individual pieces are tuned relative to these values

%     if decreasing_tension
%         dec_tension_coeff_base      = 4.6 * (3/2);  
%         valve.c_dec_tension_chordae = 1.0;     
%     else 
%         valve.dec_tension  = 0.0; 
%     end 


    % tension coefficients structure 

    % pressure / tension coefficient ratio
    % this tension coefficient is the maximum tension that a fiber can support
    tension_coeffs.pressure_tension_ratio = 0.04; 
    
    tension_coeffs.dec_tension_coeff_base = 4.6 * (3/2); 


    % tension coefficients 
    tension_coeffs.alpha_anterior             = 0.9;  % circumferential 
    tension_coeffs.beta_anterior              = 1.1;  % radial
    tension_coeffs.alpha_posterior            = 0.8;  % circumferential 
    tension_coeffs.beta_posterior             = 0.8;  % radial
    tension_coeffs.alpha_commissure           = 1.0;  % circumferential 
    tension_coeffs.beta_commissure            = 0.6;  % radial
    tension_coeffs.alpha_hoops                = 1.0;  % circumferential hoops 
    tension_coeffs.alpha_edge_connector       = 1.1;  % circumferential free edge connector 
    tension_coeffs.beta_edge_connector        = 0.01;  % circumferential free edge connector 

    % decreasing tension coefficients 
    tension_coeffs.c_circ_dec_anterior        = 1.0;  % circumferential 
    tension_coeffs.c_rad_dec_anterior         = 1.25; % radial
    tension_coeffs.c_circ_dec_posterior       = 1.0;  % circumferential 
    tension_coeffs.c_rad_dec_posterior        = 0.5;  % radial
    tension_coeffs.c_circ_dec_commissure      = 2.0;  % circumferential 
    tension_coeffs.c_rad_dec_commissure       = 1.0;  % radial 
    
    tension_coeffs.c_circ_dec_hoops           = 1.25;  % circumferential hoops
    tension_coeffs.c_rad_dec_hoops_anterior   = 0.5;  % radial hoops, anterior part 
    tension_coeffs.c_rad_dec_hoops_posterior  = 0.5;  % radial hoops, posterior part 
    tension_coeffs.c_rad_dec_hoops_commissure = 0.1;  % radial hoops, commissure part
    
    tension_coeffs.c_circ_dec_edge_connector  = 1.0;  % circumferential hoops
    tension_coeffs.c_rad_dec_edge_connector   = 1.0;  % circumferential hoops
    
    % places this many periodic rings above leaflets 
    n_rings_periodic = max(1,N/64); 
    
    % places circumferential fibers this many below hoops 
    % if the location is not already covered by leaflet 
    n_edge_connectors = max(1,N/32); 


    % Explicit commissural leaflet here 
    N_anterior = N/2; 
    angles.anterior = 5*pi/6 + pi/12; 

    % Posterior takes whatever is left 
    N_posterior = N/2;
    angles.posterior = 3*pi/6;

    N_commissure = N/4; 
    
    N_orig = N; 

    % we have added two commissural leaflets that take one fourth the total N 
    % changing the total N
    % this is a strane hack but I'm rolling with it 
    N = (3/2) * N; 

    N_per_direction   = [N_anterior/2, N_anterior/2, ... % N_anterior/4, N_anterior/4, N_anterior/4, N_anterior/4, ...  % 
                         N_commissure/2, N_commissure/2, ... 
                         N_posterior/2, N_posterior/2, ... 
                         N_commissure/2, N_commissure/2]; 

    % N_posterior/2, N_posterior/2, ... 
    %                         N_posterior/4, N_posterior/4, N_posterior/4, N_posterior/4, ... 
    %7*N_posterior/16, N_posterior/16, N_posterior/16, 7*N_posterior/16, ... % little flat center on posterior                 
    
    % store these 
    valve.N_anterior   = N_anterior; 
    valve.N_posterior  = N_posterior;
    valve.commissural_leaflets = true; 
    valve.N_commissure = N_commissure; 
    valve.N_orig       = N_orig; 
                     
                     
    % Anterior goes down then up 
    leaflet_direction = [-1, 1]; 
    % leaflet_direction = [-1, 0, 0, 1]; 
    
    % Commissure down up 
    leaflet_direction = [leaflet_direction, -1, 1]; 
    
    % Posterior goes down then up 
    leaflet_direction = [leaflet_direction, -1, 1]; 
    %leaflet_direction = [leaflet_direction, -1, 0, 0, 1]; 
    
    % Commissure down up 
    leaflet_direction = [leaflet_direction, -1, 1]; 
    
    % No offset, starting at commissure 
    leaflet_N_start = 0; 

    % changes the whole tree tension by this constant 
    tension_coeffs.tree_tension_multiplier = 0.81; 

    % Leaf tensions are all modified 
    tension_coeffs.leaf_tension_base = 1.68; 

    % Base total root tension 
    % The value 0.5905 works well on each tree when using separate solves and two leaflets 
    % Controls constant tension at the root of the tree 
    tension_coeffs.root_tension_base = 0.7; 


    
    % tree constants and parameters  
    
    % Here we take one tree half commissure width 
    % and plave between the commissure and each other leaflet 
    
    tree_n_start = N_commissure/2 + 1; 

    
    % this array determines the fraction of N_orig which each tree takes up 
    % this allows us to determine initial fractions of constants that go to each tree 
    frac_of_n_orig = [1/16; 1/16; 1/16; 1/16; ...   % anterior  
                      1/ 4; 1/ 4;             ...   % anterior and comm, comm and posterior       
                      1/16; 1/16; 1/16; 1/16; ...   % posterior
                      1/ 4; 1/ 4];                  % posterior and comm, comm and anterior
    
    % change these to manipulate individial tree coefficients 
    % for sanity reasons, these shuold mostly be one unless you have a good reason to change 
    % note that these are scaled by the fraction of the leaflet that they take up 
    k_0_1_coeff    = frac_of_n_orig .*    ... 
                     [1.6; 1.0; 1.0; 1.6; ...       % anterior  
                      1.0; 1.0;           ...       % anterior and comm, comm and posterior       
                      1.0; 1.0; 1.0; 1.0; ...       % posterior
                      1.0; 1.0]; 
                  
    k_root_coeff   = frac_of_n_orig .*    ... 
                     [1.6; 1.0; 1.0; 1.6; ...       % anterior  
                      0.9; 0.9;           ...       % anterior and comm, comm and posterior       
                      1.0; 1.0; 1.0; 1.0; ...       % posterior
                      0.9; 0.9];                    % posterior and comm, comm and anterior
                  
    % leaf coefficients scale, because we expect the lengths of the leaves to decrease 
    % with refinement of the mesh  
    % 
                                                                                    % chordae
    tension_coeffs.c_dec_chordae_leaf = (1/N)  * 2 * [1.0; 1.0; 1.0; 1.0; ...       % anterior  
                                                      1.0; 1.0;           ...       % anterior and comm, comm and posterior       
                                                      1.0; 1.0; 1.0; 1.0; ...       % posterior
                                                      1.0; 1.0];                    % posterior and comm, comm and anterior
                                                 
    % root constants do not scale, because the root 
    % should maintain a consistent length when mesh is changed 
    tension_coeffs.c_dec_chordae_root = (2/96) *     [1.0; 1.0; 1.0; 1.0; ...       % anterior  
                                                      1.0; 1.0;           ...       % anterior and comm, comm and posterior       
                                                      1.0; 1.0; 1.0; 1.0; ...       % posterior
                                                      1.0; 1.0];                    % posterior and comm, comm and anterior
                                                 
    % number of anterior trees on left 
    % for splitting up papillary muscle 
    n_trees_anterior_left = 2; 
    
    n_leaves = N_orig * frac_of_n_orig; 
    
    % this is the total leaf tension on the tree 
    % scaling by n_leaves occurs automatically when counting the number of trees 
    tension_coeffs.k_0_1  = k_0_1_coeff;
    
    % root constants,
    % actual constants, not scaled in any way 
    tension_coeffs.k_root = k_root_coeff;
    
    
    % number of trees connecting to each papillary muscle 
    trees_per_side = length(frac_of_n_orig)/2; 
    

    
    papillary_left_min_angle = -5*pi/4; 
    papillary_left_max_angle =    pi/4; 
    
    papillary_left  = get_papillary_coords(valve, left_papillary_idx,  trees_per_side,  papillary_left_min_angle,  papillary_left_max_angle); 
    
    % right takes reflected angles 
    papillary_right = get_papillary_coords(valve, right_papillary_idx, trees_per_side, -papillary_left_max_angle, -papillary_left_min_angle);
    
    % left anterior go first 
    papillary = papillary_left(:,(end - n_trees_anterior_left + 1):end); 
    
    % right, whole thing 
    papillary = [papillary, papillary_right]; 
    
    % left, remaining 
    papillary = [papillary, papillary_left(:, 1:(end - n_trees_anterior_left))]; 

end 


% default value for tree offset 
if ~exist('tree_n_start', 'var')
    tree_n_start = 1; 
end 


valve.r        = valve.skeleton.r; 

% scaling for target points 
% note that this does not include copies 
% and scaling for copies is handled by the output routine 

% scales for by mesh width for consistant total mesh force on ring 
valve.target_net_unscaled       = 8/valve.N; 

% does not scale since total number of points is constant 
valve.target_papillary_unscaled = 40/128; 

% viscoelastic damping coefficients for net, does not include copies 
valve.eta_net_unscaled = valve.target_net_unscaled/5000; 

% viscoelastic damping coefficients for root attachments, does not include copies  
valve.eta_papillary_unscaled = valve.target_papillary_unscaled/500; 

% Approximate Lagrangian mesh spacing at ring 
% Used for later splitting of springs 
% If any spring is placed at more than double this length an extra vertex is placed
valve.ds = 2*pi*valve.r / N; 


[leaflet valve] = initialize_leaflet_bead_slip(name,                                ... 
                                               N,                                   ...
                                               reflect_x,                           ... 
                                               angles,                              ...
                                               papillary,                           ... 
                                               n_leaves,                            ...
                                               tree_n_start,                        ...
                                               leaflet_direction,                   ...
                                               leaflet_N_start,                     ...
                                               N_per_direction,                     ...
                                               radial_and_circumferential,          ...  
                                               tension_coeffs,                      ... 
                                               p_0,                                 ... 
                                               n_rings_periodic,                    ...
                                               n_edge_connectors,                   ... 
                                               valve);  

valve.leaflets(1) = leaflet; 
    

% viscoelastic damping coefficients springs 
% eta, damping coeff here, is multiplied by the coefficient on the 
% associated spring 
% note that linear springs and collagen springs have vastly different constants 
% and these are tuned manually to make the dashpot constants equal order of magnitude
valve.eta_multiplier_linear   = 0; 
valve.eta_multiplier_collagen = 0; 


valve_plot(valve); 
pause(.1); 

disp('Done with initialize.'); 


