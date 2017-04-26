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