function [] = output_to_ibamr_format(valve)
    % 
    % Outputs the current configuration of the leaflets to IBAMR format
    % Spring constants are computed in dimensional form 
    % 
    %
    % Input: 
    %    base_name                  File base name
    %    L                          Outputs extra mesh to use a [-L,L]^3 cube
    %    ratio                      Ratio of pressure to nondimensionalized spring constant   
    %    params_posterior           Parameters for various leaflets 
    %    filter_params_posterior
    %    params_anterior
    %    p_physical                 Pressure in mmHg
    %    target_multiplier          Target spring strength is target_multiplier * p_physical/ratio
    %    refinement                 N/32. Mesh is this many times finer than original debug width 
    %    n_lagrangian_tracers       Places a 3d uniform mesh of tracers, this many (plus one) per side
    %    X_config_is_reference      Replaces the reference configuration with the current configuration
    %                               This attempts to remove initial transients entirely 
    % 
    % 
    % Output: 
    %    Files written in IBAMR format 
    % 
    
    N                         = valve.N; 
    base_name                 = valve.base_name; 
    L                         = valve.L; 
    posterior                 = valve.posterior; 
    anterior                  = valve.anterior; 
    tension_base              = valve.tension_base; 
    du                        = valve.anterior.du; 
    target_multiplier         = valve.target_multiplier; 
    n_lagrangian_tracers      = valve.n_lagrangian_tracers; 
    num_copies                = valve.num_copies; 
    collagen_springs_leaflet  = valve.collagen_springs_leaflet; 
    
    
    if collagen_springs_leaflet 
        error('collagen curve springs not implemented in new units'); 
    end   

    params.vertex = fopen(strcat(base_name, '.vertex'), 'w'); 
    params.spring = fopen(strcat(base_name, '.spring'), 'w'); 
    params.target = fopen(strcat(base_name, '.target'), 'w'); 
    params.inst   = fopen(strcat(base_name, '.inst'  ), 'w'); 
    
    % just make this ridiculously big for now 
    % would be better to implement some resizing but that will also clutter things up 
    params.vertices_capacity = 1000 * N^2; 
    params.vertices = zeros(3,params.vertices_capacity); 

    % keep one global index through the whole thing 
    % every time a vertex is placed this is incremented 
    % and params.vertics(:,i) contains the coordinates 
    params.global_idx = 0;
    
    % also keep an index for placing after known indices 
    params.max_idx_after_reserved_indices = 0; 
    
    % just count the number of vertices and strings throughout 
    params.total_vertices = 0; 
    params.total_springs  = 0; 
    params.total_targets  = 0; 


    % Spring constant base for targets and 
    % Approximate force is tension_base multiplied by a length element 
    k_rel = tension_base * du / num_copies; 
        
    % base rate for target spring constants
    % target constant for a single point 
    % this does not scale when the mesh is changed 
    k_target = target_multiplier * k_rel; 
    
    % No general target damping for now 
    eta = 0.0; 
    
    % the valve ring is 1d, should be halfed with doubling of mesh 
    % also set damping coefficients accordingly 
    k_target_ring = k_target; %  / refinement; 
    m_ring = 0.0; 
    eta_ring = 0.0; %sqrt(m_ring * k_target_ring); 
    
    % there are four times as many, so they get multiplied by refinement squared 
    % can also just divide by refinement because not want them to get stiffer
    k_target_net = k_target; % / refinement; 
    m_net = 0.0; 
    eta_net = 0.0; % sqrt(m_net * k_target_net);

    % output the left and right papillary as the first two vertices and targets
    
    % Critical damping for given k, mass m is 2*sqrt(m*k) 
    % Set to half critical for first test
    % 
    % m_effective_papillary = 1.0 * pi * filter_params_posterior.r^2; 
    eta_papillary         = 0.0; %sqrt(k_target/2 * m_effective_papillary); 
    
    % Approximate mesh spacing 
    % L = 2.5 = radius (of square) in sup norm, half total length of domain 
    ds = 2*L / N;
    
    
    if isfield(posterior, 'reflect_x') && posterior.reflect_x
        error('reflection not implemented in linear post slip model'); 
    end 
    if isfield(anterior, 'reflect_x') && anterior.reflect_x
        error('Something strange, should not be reflecting on anterior leaflet')
    end 
    
    
    % copies, if needed, will be placed this far down 
    if num_copies > 1
        z_offset_vals = -linspace(0, ds, num_copies); 
    else 
        z_offset_vals = 0; 
    end 
    
    % check for consistency in chordae , all data structures must match 
    if ~(    all(posterior.left_papillary   == posterior.chordae.left_papillary)   ...
          && all(posterior.right_papillary  == posterior.chordae.right_papillary))
        error('Posterior chordae are inconsistent'); 
    end 
        
    if ~(    all(anterior.left_papillary   == anterior.chordae.left_papillary)   ...
          && all(anterior.right_papillary  == anterior.chordae.right_papillary))
        error('Anterior chordae are inconsistent'); 
    end 
    
    
    % ugh this is terrible fix it it makes my head hurt by I'm tired 
    global z_offset
    
    
    for z_offset = z_offset_vals
       
        assign_papillary = true; 
         
        [params anterior] = assign_indices_vertex_target(params, anterior, assign_papillary, k_target, eta); 
        
        if ~valve.split_papillary
            assign_papillary = false; 
            posterior.left_papillary_idx  = anterior.left_papillary_idx; 
            posterior.right_papillary_idx = anterior.right_papillary_idx; 
        end 
        
        [params posterior] = assign_indices_vertex_target(params, posterior, assign_papillary, k_target, eta); 
        
        % write springs 
        params = add_springs(params, anterior,  num_copies, ds, collagen_springs_leaflet); 
        params = add_springs(params, posterior, num_copies, ds, collagen_springs_leaflet); 

        % flat part of mesh 
        r = valve.r; 
        N_ring = 2 * N;
        h = 0.0; % ring always set at zero 
        ref_frac_net = 1.0; 

        % no radial fibers, instead geodesics from the leaflet 
        radial_fibers = false; 

        % turn the polar net off for now      
        params = place_net(params, r, h, L, N_ring, radial_fibers, k_rel, k_target_net, ref_frac_net, eta_net); 
        
        % approximate geodesic continutations of fibers 
        params = place_rays(params, anterior,  ds, valve.r, L, k_rel, k_target_net, ref_frac_net, eta_net);                                       
        params = place_rays(params, posterior, ds, valve.r, L, k_rel, k_target_net, ref_frac_net, eta_net);                    


        % flat part of mesh with Cartesian coordinates
        % inner radius, stop mesh here 
        r_cartesian = r + 4*ds; 
        params = place_cartesian_net(params, r_cartesian, h, L, ds, k_rel, k_target_net, ref_frac_net, eta_net); 
 
    end 
                                            
    if n_lagrangian_tracers > 0
        double_z = false; 
        [params, total_lagrangian_placed] = place_lagrangian_tracers(params, n_lagrangian_tracers, L, double_z); 
        particles = fopen(strcat(base_name, '.particles'), 'w'); 
        fprintf(particles, '%d\n', total_lagrangian_placed); 
    end 

    % finally, write all vertices 
    params = write_all_vertices(params); 

    % and clean up files with totals 
    fclose(params.vertex); 
    fclose(params.spring); 
    fclose(params.target); 
    fclose(params.inst  ); 

    prepend_line_with_int(strcat(base_name, '.vertex'), params.total_vertices); 
    prepend_line_with_int(strcat(base_name, '.spring'), params.total_springs); 
    prepend_line_with_int(strcat(base_name, '.target'), params.total_targets); 

end 



% nest this function so it can access the z increment
% this is bad practice and should be removed 
function params = vertex_string(params, coords)
    % prints formatted string for current vertex to vertex file   
    
    % FIXME!!! 
    global z_offset
    
    fprintf(params.vertex, '%.14f\t %.14f\t %.14f\n', coords(1), coords(2), coords(3) + z_offset); 
    params.total_vertices = params.total_vertices + 1; 
end

function params = write_all_vertices(params)
    % writes all vertices to file 

    max_idx = params.global_idx; 

    for i=1:max_idx
        params = vertex_string(params, params.vertices(:,i)); 
    end 
end 


function params = spring_string(params, idx, nbr, kappa, rest_len, function_idx)
    % prints a spring format string to string file 
    if nbr <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    fprintf(params.spring, '%d\t %d\t %.14f\t %.14f', idx, nbr, kappa, rest_len); 
    
    % index for custom spring functions 
    if exist('function_idx', 'var') 
        fprintf(spring, '\t%d', function_idx); 
    end 
    
    fprintf(params.spring, '\n'); 
    
    params.total_springs = params.total_springs + 1; 
end 


function params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, num_copies)
    % 
    % Add one or more springs 
    % If the rest length is more than 2 times the specificed mesh width
    % Then it is split into multiple springs 
    % 
    % Note that the constant passed in must be the relative spring constant 
    % 

    if nbr_idx <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    N_springs = floor(rest_len / ds); 

    max_strain = .01; 
    
    X     = params.vertices(:,idx + 1); 
    X_nbr = params.vertices(:,nbr_idx + 1); 
    
    strain = (norm(X_nbr - X) - rest_len) / rest_len; 
    
    fprintf('strain = %e, idx = %d, nbr = %d\n', strain, idx, nbr_idx)
    
    if strain > max_strain 
        warning(sprintf('strain = %e, idx = %d, nbr = %d\n', strain, idx, nbr_idx)); 
    end 
    
    % Just one spring placed here 
    if N_springs <= 1 
        k_abs = k_rel / (rest_len * num_copies); 
        params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
    else 
        
        % increment of new points 
        step = (X_nbr - X)/N_springs; 
        
        X_vertices_new = zeros(3,N_springs + 1); 
        
        % set new coordinates 
        for i=0:N_springs
            X_vertices_new(:,i+1) = X + i*step; 
        end 
        
        if X ~= X_vertices_new(:,1)
            error('First vertex must be equal to X'); 
        end 
        
        if X_nbr ~= X_vertices_new(:,N_springs + 1)
            error('Last vertex must be equal to X_nbr'); 
        end
        
        
        % first index is always first index suppied 
        idx_tmp = idx; 
        
        for i=1:N_springs 
            
            % place the upper new vertex if it is not placed already 
            if i<N_springs 
                params.vertices(:,params.global_idx + 1) = X_vertices_new(:,i+1); 
                
                % new vertex's index is current (zero indexed) global index 
                nbr_idx_tmp = params.global_idx; 

                params.global_idx = params.global_idx + 1;
                
            else 
                % last spring upper limit is the previous neighbor 
                nbr_idx_tmp = nbr_idx; 
            end 
            
            % global indices placed in order by convention 
            min_idx = min(idx_tmp, nbr_idx_tmp); 
            max_idx = max(idx_tmp, nbr_idx_tmp); 
        
            % Current length 
            L = norm(X_vertices_new(:,i+1) - X_vertices_new(:,i)); 
            
            % Rest length determined by strain 
            R = L / (strain + 1); 
            
            % Absolute spring constant must be used in spring file 
            k_abs = k_rel / (R * num_copies); 
            
            % Finally, write the spring string 
            params = spring_string(params, min_idx, max_idx, k_abs, R); 
            
            % lower index is always previous upper index 
            idx_tmp = nbr_idx_tmp; 
            
        end 
    
    end 

end 


function params = target_string(params, idx, kappa, eta)
    % prints a target format string to target file 
    
    if exist('eta', 'var') && (eta > 0.0)
        fprintf(params.target, '%d\t %.14f\t %.14f\n', idx, kappa, eta);
    else
        fprintf(params.target, '%d\t %.14f\n', idx, kappa);
    end 
    params.total_targets = params.total_targets + 1; 
end 


function [] = prepend_line_with_int(file_name, val)
    % Adds a single line to the file with given name
    % at the beginning with the integer val 
    % writes a temp file then calls cat 

    
    write = fopen('temp.txt', 'w'); 
    fprintf(write, '%d\n', val); 
    fclose(write); 
    
    file_temp = strcat(file_name, '.tmp'); 
    
    system(sprintf('cat temp.txt %s > %s', file_name, file_temp));
    movefile(file_temp, file_name); 
    system('rm temp.txt'); 

end 


function [params leaflet] = assign_indices_vertex_target(params, leaflet, assign_papillary, k_target, eta)
    % 
    % Assigns global indices to the leaflet and chordae 
    % Includes all boundary condition points and internal 
    % 
    % 
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    
    % Unpack needed data 
    X                 = leaflet.X; 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    is_bc             = leaflet.is_bc; 
    chordae           = leaflet.chordae;

    % Keep track of vector index in 1d array 
    leaflet.indices_global = nan * zeros(j_max, k_max); 
    
    if isfield(leaflet, 'chordae') && ~isempty(chordae)
        [m N_chordae] = size(leaflet.chordae.C_left); 
    else 
        error('Leaflet without chordae not implemented'); 
    end 
    
    leaflet.chordae.indices_global_left  = nan * zeros(N_chordae, 1); 
    leaflet.chordae.indices_global_right = nan * zeros(N_chordae, 1); 
    
    % assign indices and vertex coordinates for the papillary tips
    if assign_papillary
            
        params.vertices(:,params.global_idx + 1) = leaflet.left_papillary; 
        leaflet.left_papillary_idx               = params.global_idx; 
        params                                   = target_string(params, params.global_idx, k_target, eta);     
        params.global_idx                        = params.global_idx + 1; 

        params.vertices(:,params.global_idx + 1) = leaflet.right_papillary; 
        leaflet.right_papillary_idx              = params.global_idx; 
        params                                   = target_string(params, params.global_idx, k_target, eta);     
        params.global_idx                        = params.global_idx + 1;   
                
    else 
        % If they are not already written 
        % their indices must aleady be assignedd 
        if (~isfield(leaflet, 'left_papillary_idx')) || (~isfield(leaflet, 'right_papillary_idx'))
            error('must have papillary indices assigned before call, or assign them here'); 
        end 
    end 
    

    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point written to the file 
            if is_internal(j,k) || is_bc(j,k)
                
                params.vertices(:,params.global_idx + 1) = X(:,j,k); 
                leaflet.indices_global(j,k) = params.global_idx; 
                
                % if on boundary, this is a target point 
                if is_bc(j,k)
                    if exist('eta', 'var')
                        params = target_string(params, params.global_idx, k_target, eta);     
                    else
                        params = target_string(params, params.global_idx, k_target);     
                    end 
                end 
                
                params.global_idx = params.global_idx + 1;
            end 
        end 
    end
    
    for left_side = [true false];  

        for i=1:N_chordae

            if left_side 
                params.vertices(:,params.global_idx + 1) = chordae.C_left(:,i); 
                leaflet.chordae.indices_global_left(i)   = params.global_idx; 
            else 
                params.vertices(:,params.global_idx + 1) = chordae.C_right(:,i); 
                leaflet.chordae.indices_global_right(i)  = params.global_idx;                 
            end 
            
            params.global_idx = params.global_idx + 1;
        end
    end 

end 
 

function params = add_springs(params, leaflet, num_copies, ds, collagen_spring)

    params = add_leaflet_springs(params, leaflet, num_copies, ds, collagen_spring); 
    params = add_chordae_tree_springs(params, leaflet, num_copies, ds, collagen_spring); 

end 


function params = add_leaflet_springs(params, leaflet, num_copies, ds, collagen_spring)
                      
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    % Unpack needed data 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    is_bc             = leaflet.is_bc; 
    chordae           = leaflet.chordae;
    chordae_idx_left  = leaflet.chordae_idx_left; 
    chordae_idx_right = leaflet.chordae_idx_right; 

    R_u = leaflet.R_u;
    k_u = leaflet.k_u;
    R_v = leaflet.R_v;
    k_v = leaflet.k_v;

    R_free_edge_left   = leaflet.R_free_edge_left;
    k_free_edge_left   = leaflet.k_free_edge_left;
    R_free_edge_right  = leaflet.R_free_edge_right;
    k_free_edge_right  = leaflet.k_free_edge_right; 
    
    
    if collagen_spring
        error('not implemeneted'); 
        function_idx = 1;
    end 
    
    if isfield(leaflet, 'chordae') && ~isempty(chordae)
        [m N_chordae] = size(leaflet.chordae.C_left); 
    else 
        error('Leaflet without chordae not implemented'); 
    end 
    
    
    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point written to the file 
            if is_internal(j,k) || is_bc(j,k)
                
                % global index of current opint 
                idx = leaflet.indices_global(j,k); 
                
                % Connect to the left papillary or chordae tree 
                if chordae_idx_left(j,k)
                    
                    % current node has a chordae connection
                    
                    % chordae_idx tells what free edge index to use 
                    i = chordae_idx_left(j,k); 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_left(j,k) + N_chordae; 

                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  
                    
                    nbr_idx = chordae.indices_global_left(idx_chordae); 

                    rest_len = R_free_edge_left(i); 

                    k_rel = k_free_edge_left(i); 
                    
                    if collagen_spring 
                        % relative constant here
                        % other parameters coded into function 
                        kappa = k_rel * chordae.k_0; 
                        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx); 
                    else 
                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, num_copies); 
                    end 

                end 
                
                % Connect to the right papillary or chordae tree 
                if chordae_idx_right(j,k)
                    
                    % current node has a chordae connection
                    
                    i = chordae_idx_right(j,k); 
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_right(j,k) + N_chordae; 

                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  

                    nbr_idx = chordae.indices_global_right(idx_chordae); 
                        
                    rest_len = R_free_edge_right(i); 

                    k_rel = k_free_edge_right(i); 

                    if collagen_spring 
                        % relative constant here
                        % other parameters coded into function 
                        kappa = k_rel * chordae.k_0; 
                        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx); 
                    else
                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, num_copies);
                    end 

                end 
                                
                
                % springs in leaflet, only go in up direction 
                j_nbr = j + 1; 
                k_nbr = k; 
                if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))
                    
                    % since always moving in up direction, j_spr = j, k_spr = k
                    rest_len = R_u(j,k); 
                    k_rel    = k_u(j,k); 
                    
                    nbr_idx = leaflet.indices_global(j_nbr,k_nbr);
                    
                    if collagen_spring
                        kappa = alhpa * k_rel;         
                        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx); 
                    else 
%                         k_abs = k_rel / (rest_len * num_copies); 
%                         params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
                        
                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, num_copies);
                    end 

                end 
                
                % springs in leaflet, only go in up direction 
                j_nbr = j; 
                k_nbr = k + 1; 
                if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))
                    
                    % since always moving in up direction, j_spr = j, k_spr = k
                    rest_len = R_v(j,k); 
                    k_rel    = k_v(j,k); 
                    
                    nbr_idx = leaflet.indices_global(j_nbr,k_nbr);
                    
                    if collagen_spring
                        kappa = beta * k_rel;         
                        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx); 
                    else 
                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, num_copies);
                    end 

                end 
                
            end 
        end 
    end
end 



function params = add_chordae_tree_springs(params, leaflet, num_copies, ds, collagen_spring)

    % Adds chordae tree to IBAMR format files 
    % No targets here, so files and and count not included 
                        
    if ~isfield(leaflet, 'chordae') || isempty(leaflet.chordae)
        error('cannot place chordae with empty array'); 
    end 
    
    chordae         = leaflet.chordae; 

    [m N_chordae] = size(chordae.C_left); 
    
    if collagen_spring
        error('not implemeneted'); 
        function_idx = 1;
    end
    
    
    for left_side = [true false];  
        
        if left_side
            indices_global = chordae.indices_global_left; 
        else 
            indices_global = chordae.indices_global_right; 
        end 
        
        for i=1:N_chordae

            idx = indices_global(i); 
            
            % place the spring which goes to the parent 
            parent = floor(i/2); 
            
            if parent == 0
                
                % zero index means papillary muscles
                if left_side 
                    nbr_idx = leaflet.left_papillary_idx; 
                else 
                    nbr_idx = leaflet.right_papillary_idx; 
                end
                
            else
                nbr_idx = indices_global(parent); 
            end 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr rest_len k_rel] = get_nbr_chordae(leaflet, i, parent, left_side); 
            
            if collagen_spring 
                kappa = k_rel * k_val; 
                % list nbr index first because nbr is parent and has lower index
                params = spring_string(params, nbr_idx, idx, kappa, rest_len, function_idx);            
            else 
                % list nbr index first because nbr is parent and has lower index
                params = place_spring_and_split(params, nbr_idx, idx, k_rel, rest_len, ds, num_copies);
                
            end 
        end 
        
    end 

end 
                        
                        
function params = place_net(params, r, h, L, N, radial_fibers, k_rel, k_target, ref_frac, eta)
    % 
    % Places a polar coordinate mesh in a box 
    % Starts with N points evenly spaced at radius R
    % Spacing ds is determined by this value 
    % Evenly spaced in R, theta spacing increases as radius increases 
    %
    % Mesh is placed with a rectangular topology
    % Springs everywhere internal 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     N                Number of points placed on the initial circle (valve ring)
    %     spring           Files for writing, must be open already 
    %     vertex 
    %     target
    %     inst 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    % mesh spacing on the valve ring 
    ds = 2*pi / N; 

    % maximum number of points in radial direction 
    % this gets a sqrt(2) because we are placing the net in a square 
    if radial_fibers
        % place rings to edge 
        M = floor(sqrt(2) * L / ds); 
    else 
        % only include those full rings which fit in the domain 
        M = floor( (L-r) / ds); 
    end 
        
    points = zeros(3,N,M);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(N,M); 
    
    instrument_idx = 0; 
    
    rad = r; 
    
    % compute vertex positions and add to array 
    for k=1:M
        for j=1:N
            
            theta = (j-1) * ds; 
            coords_horiz = [rad*cos(theta); rad*sin(theta)]; 
            
            % if one norm is less than L, then the point is within the domain  
            if norm(coords_horiz, inf) < L 
                points(:,j,k) = [coords_horiz; h]; 
                indices_global(j,k) = params.global_idx; 
                params.vertices(:,params.global_idx + 1) = points(:,j,k); 
                
                % every valid vertex is a target point here 
                if exist('eta', 'var')
                    params = target_string(params, params.global_idx, k_target, eta);     
                else
                    params = target_string(params, params.global_idx, k_target);     
                end 

                params.global_idx = params.global_idx + 1;                   

            else 
                points(:,j,k) = NaN * ones(3,1);
                indices_global(j,k) = NaN; 
            end     
            
        end 
        
        rad = rad + ds; 
    end 

    
    % write the instrument file header here 
    fprintf(params.inst, '1   # num meters in file\n'); 
    fprintf(params.inst, 'meter_0   # name\n'); 
    fprintf(params.inst, '%d  # number of meter points\n', N); 
    
    
    % below the first possible point 
    idx = -1; 
    
    for k=1:M
        for j=1:N
            
            % just ignore the nan 
            if ~isnan(indices_global(j,k))
 
                last_idx = idx; 
                idx = indices_global(j,k); 
                
                % instrument file on 
                if k == 1
                    fprintf(params.inst, '%d \t0 \t %d\n', idx, instrument_idx); 
                    instrument_idx = instrument_idx + 1; 
                end 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 

                % check up directions for springs 
                if (j+1) < N
                    if ~isnan(indices_global(j+1,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j+1,k)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(j+1,k); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
                    end 
                end 
                
                % don't forget the periodic direction in j
                if (j+1) == N
                   % need to make sure that the 1,k point is also not a NaN  
                   if ~isnan(indices_global(1,k)) 
                       rest_len = ref_frac * norm(points(:,j,k) - points(:,1,k)); 
                       k_abs = k_rel / rest_len;
                       nbr_idx = indices_global(1,k); 
                       params = spring_string(params, nbr_idx, idx, k_abs, rest_len); 
                   end 
                end 
                
                % check up directions for springs 
                if radial_fibers
                    if (k+1) < M
                        if ~isnan(indices_global(j,k+1))
                            rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                            k_abs = k_rel / rest_len; 
                            nbr_idx = indices_global(j,k+1); 
                            params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
                        end 
                    end 
                end 
                
                % no periodic direction in radial direction (k)
                
            end 
           
        end 
    end  

end 


function params = place_rays(params, leaflet, ds, r, L, k_rel, k_target, ref_frac, eta)
    % 
    % Places rays of fibers emenating from the leaflet 
    % Angle of rays makes them (roughly) geodesics 
    % on the valve surface and plane surface  
    % 
    % Spacing is determined by reflection 
    %
    % Mesh is placed with a rectangular topology
    % Springs everywhere internal 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     N                Number of points placed on the initial circle (valve ring)
    %     spring           Files for writing, must be open already 
    %     vertex 
    %     target
    %     inst 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    X           = leaflet.X; 
    j_max       = leaflet.j_max; 
    k_max       = leaflet.k_max;
    is_bc       = leaflet.is_bc; 
    is_internal = leaflet.is_internal; 
    h           = 0.0;        % always place at origin 
    
    
    if ~isfield(leaflet, 'indices_global')
        error('Must place leaflets before placing rays'); 
    end 
    
    
    for j = 1:j_max
        for k = 1:k_max
            if is_bc(j,k)
            
                pt_ring = X(:,j,k); 

                % only get a fiber if the previous point is included in the leaflet  
                neighbors = []; 
                j_nbr = j-1; 
                k_nbr = k; 
                if (j_nbr > 0) &&  (k_nbr > 0) && is_internal(j_nbr,k_nbr)
                    neighbors = [X(:,j-1,k), neighbors] ; 
                end
                
                j_nbr = j; 
                k_nbr = k-1; 
                if (j_nbr > 0) &&  (k_nbr > 0) && is_internal(j_nbr,k_nbr)
                    neighbors = [X(:,j,k-1), neighbors] ; 
                end        

                for x = neighbors 

                    % find the initial reflected point 
                    val = get_geodesic_continued_point(x, pt_ring, r, h); 

                    % each point moves by this much from the initial point 
                    increment = val - x; 

                    % just zero this component, they are both near 3 but maybe not exactly 
                    increment(3) = 0.0; 
                    
                    % set to proper mesh width regardless of other spacing 
                    increment = ds * increment / norm(increment); 


                    % point = val; 
                    point_prev = pt_ring; 
                    point = pt_ring + increment; 
                    nbr_idx = leaflet.indices_global(j,k); 

                    % just keep adding until points leave the domain 
                    while norm(point(1:2), inf) < L   

                        % grab the index 
                        idx = params.global_idx;

                        % place point 
                        params.vertices(:,params.global_idx + 1) = point; 
                        
                        % it's a target too 
                        if exist('eta', 'var')
                            params = target_string(params, idx, k_target, eta);     
                        else
                            params = target_string(params, idx, k_target);     
                        end 

                        rest_len = ref_frac * norm(point - point_prev); 
                        k_abs = k_rel / rest_len;

                        params = spring_string(params, nbr_idx, idx, k_abs, rest_len);

                        point_prev = point; 
                        point      = point + increment; 
                        params.global_idx = params.global_idx + 1; 
                        nbr_idx    = idx; 
                    end 

                end 
            end 
        end 
    end 
end 


function [val] = get_geodesic_continued_point(x, pt_ring, r, h)
    %
    % Takes a point inside the the valve ring
    % Returns a point which allows a geodesic continuation 
    % Taking a segment from the ring point to the wall creates a geodesic 
    % 
    % Input: 
    %     x           Coordinates 
    %     pt_ring     Point on the valve ring to reflect relative to 
    %     r           Radius of valve ring 
    %     h           height of valve ring 
    %

    tol = 1e5 * eps; 
    
    if abs(pt_ring(3) - h) > tol 
        error('Initial ring point is not at the right height'); 
    end 
    
    if abs(norm(pt_ring(1:2)) - r) > tol 
        error('Initial ring point is not at the correct radius'); 
    end 
    
    
    % angle needed to send the ring point to the x axis 
    theta = atan2(pt_ring(2), pt_ring(1)); 
    
    % rotate the point of inter
    val = rotation_matrix_z(-theta) * x; 
    
    % ring point to origin, val to near origin to be reflected 
    val = val - [r; 0; h]; 

    % rotate val into the z = 0 plane, inside the transformed ring 
    phi = atan2(val(3), val(1)); 
    val = rotation_matrix_y( -(phi - pi)) * val;
    
    if abs(val(3)) > tol
        error('should have landed in z=0 plane here...');
    end 

    % reflect, this would be the geodesic point if the system was flat 
    val = -val; 
    
    % send back to original coordinates 
    % val = rotation_matrix_y(phi - pi) * val;
    val = val + [r; 0; h]; 
    val = rotation_matrix_z(theta) * val; 
    
    if abs(val(3) - h) > tol 
        error('rotated value is not in plane'); 
    end 
    
    if norm(val(1:2)) <= r 
        error('rotated value must be out of the plane'); 
    end 
    
end 






function params = place_cartesian_net(params, r, h, L, ds, k_rel, k_target, ref_frac, eta)
    % 
    % Places a cartesian coordinate mesh in a box 
    % This is to avoid issues with the polar mesh at the edge 
    % Mesh is placed on [-L + ds/2, L-ds/2]
    % 
    % Springs everywhere internal in 2d arrangement 
    % All points are made targets
    % 
    % Input
    %     r                Internal radius, it is advisible to make this larger than the valve ring 
    %     h                Height of valve ring 
    %     L                Box width -- point must have max norm less than L to be included 
    %     ds               Mesh is placed starting ds from boundary 
    %     params.global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    
    N = ceil(2*L / ds);  
    points = zeros(3,N,N);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(N,N);     
    
    % This loop should be one indexed, 
    % because we want to start not at the edge but ds in 
    for k=1:N
        for j=1:N
            
            coords_horiz = [ (j-1)*ds - L + ds/2; (k-1)*ds - L + ds/2]; 
            
            % if one norm is less than L, then the point is within the domain  
            if (norm(coords_horiz, inf) < L) && (norm(coords_horiz,2) > r) 
                
                points(:,j,k) = [coords_horiz; h]; 
                indices_global(j,k) = params.global_idx; 
                params.vertices(:,params.global_idx + 1) = points(:,j,k);

                % every valid vertex is a target point here 
                if exist('eta', 'var')
                    params = target_string(params, params.global_idx, k_target, eta);     
                else
                    params = target_string(params, params.global_idx, k_target);     
                end 

                params.global_idx = params.global_idx + 1; 
                
                
            else 
                points(:,j,k) = NaN * ones(3,1);
                indices_global(j,k) = NaN; 
            end     
            
        end 
        
    end 

    
    % below the first possible point 
    idx = -1; 
    
    for k=1:N
        for j=1:N
            
            % just ignore the NaNs 
            if ~isnan(indices_global(j,k))
 
                last_idx = idx; 
                idx = indices_global(j,k); 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 

                % check up directions for springs 
                if (j+1) <= N
                    if ~isnan(indices_global(j+1,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j+1,k)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(j+1,k); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
                    end 
                end 
                                
                % check up directions for springs 
                if (k+1) <= N
                    if ~isnan(indices_global(j,k+1))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                        k_abs = k_rel / rest_len; 
                        nbr_idx = indices_global(j,k+1); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len); 
                    end 
                end 
                
            end 
           
        end 
    end 

end 



function [params, total_lagrangian_placed] = place_lagrangian_tracers(params, n_lagrangian_tracers, L, double_z)
    % Places a uniform cartesian mesh of lagrangian particle tracers 
    % Simple lopp implementation 
    %
    %     params.global_idx               Running totals  
    %     total_vertices 
    %     vertex                   vertex file for writing  
    %     n_lagrangian_tracers
    %     L                        mesh placed in L/2

    % junk hack exit 
    if n_lagrangian_tracers == 0
        total_lagrangian_placed = 0; 
        return; 
    end 
    
    dx = L / n_lagrangian_tracers; 
    total_lagrangian_placed = 0; 
    
    z_extra_offset = 0.0; 
    
    n_z_dir = n_lagrangian_tracers; 
    z_min = -L/2; 
    if double_z
        n_z_dir = 2*n_z_dir; 
        z_min = -L; 
    end 
    
    for i = 0:n_lagrangian_tracers
        for j = 0:n_lagrangian_tracers
            for k = 0:n_z_dir; 
                
                x = i * dx - L/2; 
                y = j * dx - L/2; 
                z = k * dx + z_min + z_extra_offset; 
                
                params.vertices(:,params.global_idx + 1) = [x y z]; 
    
                total_lagrangian_placed = total_lagrangian_placed + 1; 
                params.global_idx = params.global_idx + 1; 
            end 
        end 
    end 

end 





