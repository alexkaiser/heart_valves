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
    
    N                           = valve.N; 
    base_name                   = valve.base_name; 
    L                           = valve.L; 
    tension_base                = valve.tension_base; 
    target_net                  = valve.target_net; 
    target_papillary            = valve.target_papillary; 
    n_lagrangian_tracers        = valve.n_lagrangian_tracers; 
    collagen_constitutive       = valve.collagen_constitutive; 

    % if this is true, all partition gets ignored except right at valve
    % ring 
    if ~isfield(valve, 'in_heart')
        in_heart = false; 
    else 
        in_heart = valve.in_heart; 
    end 
    
    if in_heart 
        n_lagrangian_tracers = 0; 
    end 
    
    if ~valve.split_papillary
        error('Must have split papillary locations in current implementation.'); 
    end

    params.vertex    = fopen(strcat(base_name, '.vertex'), 'w'); 
    params.spring    = fopen(strcat(base_name, '.spring'), 'w'); 
    params.target    = fopen(strcat(base_name, '.target'), 'w'); 
    params.inst      = fopen(strcat(base_name, '.inst'), 'w'); 
    params.papillary = fopen(strcat(base_name, '.papillary'), 'w'); 
    
    % just make this ridiculously big for now 
    % would be better to implement some resizing but that will also clutter things up 
    params.vertices_capacity = 100 * N^2; 
    params.vertices = zeros(3,params.vertices_capacity); 

    % keep one global index through the whole thing 
    % every time a vertex is placed this is incremented 
    % and params.vertics(:,i) contains the coordinates 
    params.global_idx = 0;
    
    % also keep an index for placing after known indices 
    params.max_idx_after_reserved_indices = 0; 
    
    % just count the number of vertices and strings throughout 
    params.total_vertices  = 0; 
    params.total_springs   = 0; 
    params.total_targets   = 0; 
    params.total_papillary = 0; 

    % keep a single parameter for outputting copies 
    params.z_offset = 0; 

    % box sizes 
    params.x_min = -L; 
    params.x_max =  L; 
    params.y_min = -L; 
    params.y_max =  L; 
    params.z_min =  3.0 - 4*L; 
    params.z_max =  3.0; 
    
    if isfield(valve.skeleton, 'ring_center')
        ring_center = valve.skeleton.ring_center
        params.x_min = params.x_min + valve.skeleton.ring_center(1); 
        params.x_max = params.x_max + valve.skeleton.ring_center(1); 
        params.y_min = params.y_min + valve.skeleton.ring_center(2); 
        params.y_max = params.y_max + valve.skeleton.ring_center(2); 
        params.z_min = params.z_min + valve.skeleton.ring_center(3); 
        params.z_max = params.z_max + valve.skeleton.ring_center(3); 
    else
        valve.skeleton.ring_center = zeros(3,1); 
    end 
    params.ring_center = valve.skeleton.ring_center; 
    
    % parameters for scaling of other constants 
    params.num_copies = valve.num_copies; 
    params.eta_multiplier_linear   = valve.eta_multiplier_linear; 
    params.eta_multiplier_collagen = valve.eta_multiplier_collagen; 
    
    % parameters for output flags 
    params.output = valve.output; 
    
    % Spring constant base for targets and 
    % Approximate force is tension_base multiplied by a length element 
    du = 1/N; 
    k_rel = tension_base * du; 
        

    % papillary target constants 
    % this does not scale when the mesh is changed 
    k_target_papillary = target_papillary; 
    eta_papillary      = valve.eta_papillary; 
    
    k_target_net       = target_net; 
    eta_net            = valve.eta_net; 

    
    % Lagrangian mesh spacing 
    ds = valve.ds;
    
    
    if isfield(valve, 'kappa_cross_layer_multipler') && (valve.kappa_cross_layer_multipler ~= 0) && (params.num_copies > 1)
        params.cross_layer_on          = true; 
        params.kappa_cross_layer       = valve.kappa_cross_layer_multipler * tension_base; 
        params.rest_len_cross_layer    = ds; 
        params.total_per_layer         = nan; 
        params.min_idx_for_cross_layer = nan; 
        params.max_idx_for_cross_layer = nan; 
    else 
        params.cross_layer_on          = false; 
    end 
    
    
    
    % copies, if needed, will be placed this far down 
    if params.num_copies > 1
        z_offset_vals = -linspace(0, ds, params.num_copies); 
    else 
        z_offset_vals = 0; 
    end  
    
    % print the increment for systolic motion
    % this is the negative of the motion that occurred 
    % from the systolic to diastolic before 
    diastolic_increment = valve.diastolic_increment; 
    fprintf(params.papillary, '%.14f\t %.14f\t %.14f\n', diastolic_increment(1), diastolic_increment(2), diastolic_increment(3)); 
    times = valve.papillary_movement_times; 
    if length(times) ~= 5
        error('Must have five times for movement'); 
    end 
    fprintf(params.papillary, '%.14f\t %.14f\t %.14f %.14f %.14f\n', times(1), times(2), times(3), times(4), times(5)); 
    
    for copy = 1:params.num_copies
        
        params.copy = copy; 
        
        
        if params.cross_layer_on
            params.min_idx_for_cross_layer = params.global_idx; 
        end 
        
        params.z_offset = z_offset_vals(copy); 
        first_idx = params.global_idx + 1; 
        
        for i=1:length(valve.leaflets)
            j_max = valve.leaflets(i).j_max; 
            k_max = valve.leaflets(i).k_max; 
            valve.leaflets(i).indices_global = nan * zeros(j_max, k_max);  
        
            % allocate chordae indexing arrays here 
            % to stop dissimilar struct assignment error 
            for tree_idx = 1:valve.leaflets(i).num_trees
                C = valve.leaflets(i).chordae(tree_idx).C; 
                [m N_chordae] = size(C);
                valve.leaflets(i).chordae(tree_idx).idx_root = nan;                 
                valve.leaflets(i).chordae(tree_idx).indices_global = nan * zeros(N_chordae, 1);
            end 
        end 
        
        for i=1:length(valve.leaflets)
            [params valve.leaflets(i)] = assign_indices_vertex_target(params, valve.leaflets(i), k_target_net, k_target_papillary, eta_net, eta_papillary); 
        end 
        
        for i=1:length(valve.leaflets)
            params = add_springs(params, valve.leaflets(i), ds, collagen_constitutive); 
        end 

        if params.cross_layer_on
            params.max_idx_for_cross_layer = params.global_idx; 
            if copy > 1 
                params = place_cross_layer_springs(params); 
            end 
        end 
        
        % flat part of mesh 
        r = valve.r; 
        ref_frac_net = 1.0;

        if ~in_heart
            for i=1:length(valve.leaflets)
                if length(valve.leaflets) ~= 1 
                    error('Only one leaflet version currently supported'); 
                end 
                params = place_net(params, valve.leaflets(i), ds, r, L, k_rel, k_target_net, ref_frac_net, eta_net); 
            end 

            % approximate geodesic continutations of fibers 
            for i=1:length(valve.leaflets)
                params = place_rays(params, valve.leaflets(i), ds, L, k_rel, k_target_net, ref_frac_net, eta_net);
            end 

            % flat part of mesh with Cartesian coordinates
            % inner radius, stop mesh here 
            r_extra = 4*ds; 
            for i=1:length(valve.leaflets)
                if length(valve.leaflets) ~= 1 
                    error('Only one leaflet version currently supported'); 
                end 
                params = place_cartesian_net(params, valve.leaflets(i), r_extra, L, ds, k_rel, k_target_net, ref_frac_net, eta_net); 
            end 
            
        else 
            
            % pass L=r to get only one ring placed 
            params = place_net(params, valve.leaflets(i), ds, r, r, k_rel, k_target_net, ref_frac_net, eta_net); 
            
        end 
        
        
        
        % first time through, count all included indices 
        if params.cross_layer_on && (copy == 1)
            params.total_per_layer = params.global_idx; 
        end 
        
        % adjust for offset 
        last_idx = params.global_idx + 1; 
        params.vertices(3,first_idx:last_idx) = params.vertices(3,first_idx:last_idx) + params.z_offset; 
        
    end 

    if n_lagrangian_tracers > 0
        double_z = true; 
        [params, total_lagrangian_placed] = place_lagrangian_tracers(params, n_lagrangian_tracers, double_z); 
        particles = fopen(strcat(base_name, '.particles'), 'w'); 
        fprintf(particles, '%d\n', total_lagrangian_placed); 
        fclose(particles); 
    end 
    
    % finally, write all vertices 
    params = write_all_vertices(params); 

    % and clean up files with totals 
    fclose(params.vertex   ); 
    fclose(params.spring   ); 
    fclose(params.target   ); 
    fclose(params.inst     ); 
    fclose(params.papillary); 

    prepend_line_with_int(strcat(base_name, '.vertex'), params.total_vertices); 
    prepend_line_with_int(strcat(base_name, '.spring'), params.total_springs); 
    prepend_line_with_int(strcat(base_name, '.target'), params.total_targets); 
    prepend_line_with_int(strcat(base_name, '.papillary'), params.total_papillary); 

end 


function params = vertex_string(params, coords)
    % prints formatted string for current vertex to vertex file   
    
    fprintf(params.vertex, '%.14f\t %.14f\t %.14f\n', coords(1), coords(2), coords(3)); 
    params.total_vertices = params.total_vertices + 1; 
end

function params = write_all_vertices(params)
    % writes all vertices to file 

    max_idx = params.global_idx; 
    
    debug = true; 
    if debug 
        min_x = min(params.vertices(1,1:max_idx)) 
        max_x = max(params.vertices(1,1:max_idx)) 
        min_y = min(params.vertices(2,1:max_idx)) 
        max_y = max(params.vertices(2,1:max_idx)) 
        min_z = min(params.vertices(3,1:max_idx)) 
        max_z = max(params.vertices(3,1:max_idx)) 
    end 
    
    for i=1:max_idx
        params = vertex_string(params, params.vertices(:,i)); 
    end 
end 


function params = spring_string(params, idx, nbr, kappa, rest_len, function_idx, output)
    % prints a spring format string to string file 
    if nbr <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    fprintf(params.spring, '%d\t %d\t %.14f\t %.14f', idx, nbr, kappa/params.num_copies, rest_len); 
    
    % index for custom spring functions 
%     if ~exist('function_idx', 'var') 
%         function_idx = 0; 
%     end 
   
    fprintf(params.spring, '\t %d', function_idx); 
    
%     if function_idx == 0
%         eta = params.eta_multiplier_linear   * kappa / params.num_copies; 
%     elseif function_idx == 1
%         eta = params.eta_multiplier_collagen * kappa / params.num_copies; 
%     else 
%         error('Only linear (default) and collagen function indices implemented'); 
%     end 
%    
%     fprintf(params.spring, '\t %.14f', eta); 


    fprintf(params.spring, '\t # %d', output); 

    fprintf(params.spring, '\n'); 

    params.total_springs = params.total_springs + 1; 
end 


function params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output)
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
    
    if collagen_spring
        function_idx = 1;
    else 
        function_idx = 0;
    end 
    
    if k_rel == 0 
        warning(sprintf('Zero strength spring on idx,nbr = %d,%d, not placed.', idx, nbr_idx)); 
        return; 
    end 

%    max_strain = .01; 
    
    X     = params.vertices(:,idx + 1); 
    X_nbr = params.vertices(:,nbr_idx + 1); 
    L     = norm(X_nbr - X); 
    
    % Can use either rest length or current length to determine number of springs 
    % N_springs = floor(rest_len / ds); 
    N_springs = floor(L / ds); 
    
    strain = (L - rest_len) / rest_len; 
    
    % fprintf('strain = %e, idx = %d, nbr = %d\n', strain, idx, nbr_idx)
    
%     if strain > max_strain 
%         warning(sprintf('strain = %e, idx = %d, nbr = %d\n', strain, idx, nbr_idx)); 
%     end 
    
    % Just one spring placed here 
    if N_springs <= 1 
        if collagen_spring
            % Scaling constant here does not change with rest lengths 
            k_col = k_rel; 

            % Finally, write the spring string 
            params = spring_string(params, idx, nbr_idx, k_col, rest_len, function_idx, output); 
        else 
            k_abs = k_rel / rest_len; 
            params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output); 
        end 
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
            
            if collagen_spring
                % Scaling constant here does not change with rest lengths 
                k_col = k_rel; 

                % Finally, write the spring string 
                params = spring_string(params, min_idx, max_idx, k_col, R, function_idx, output);
            
            else 
                % Absolute spring constant must be used in spring file 
                k_abs = k_rel / R; 

                % Finally, write the spring string 
                params = spring_string(params, min_idx, max_idx, k_abs, R, function_idx, output);
            end 
            
            % lower index is always previous upper index 
            idx_tmp = nbr_idx_tmp; 
            
        end 
    
    end 

end 


function params = target_string(params, idx, kappa, eta)
    % prints a target format string to target file 
    
    if exist('eta', 'var') && (eta > 0.0)
        fprintf(params.target, '%d\t %.14f\t %.14f\n', idx, kappa/params.num_copies, eta/params.num_copies);
    else
        fprintf(params.target, '%d\t %.14f\n', idx, kappa/params.num_copies);
    end 
    params.total_targets = params.total_targets + 1; 
end 


function params = papillary_string(params, idx, coords)
    % prints a papillary format string to target file 
    fprintf(params.papillary, '%d\t %.14f\t %.14f %.14f\n', idx, coords(1), coords(2), coords(3) + params.z_offset);
    params.total_papillary = params.total_papillary + 1; 
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


function [params leaflet] = assign_indices_vertex_target(params, leaflet, k_target_net, k_target_papillary, eta_net, eta_papillary)
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
    num_trees         = leaflet.num_trees; 
    chordae           = leaflet.chordae; 

    % Keep track of vector index in 1d array 
   

    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point written to the file 
            if is_internal(j,k) || is_bc(j,k)
                
                params.vertices(:,params.global_idx + 1) = X(:,j,k); 
                leaflet.indices_global(j,k) = params.global_idx; 
                
                % if on boundary, this is a target point 
                if is_bc(j,k)
                    if exist('eta_net', 'var')
                        params = target_string(params, params.global_idx, k_target_net, eta_net);     
                    else
                        params = target_string(params, params.global_idx, k_target_net);     
                    end 
                end 
                
                params.global_idx = params.global_idx + 1;
            end 
        end 
    end
    
    
    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        [m N_chordae] = size(C);         

        % always place root index first 
        params.vertices(:,params.global_idx + 1) = leaflet.chordae(tree_idx).root; 
        leaflet.chordae(tree_idx).idx_root = params.global_idx;                 
        
        % root is always a boundary condition 
        if exist('eta_papillary', 'var')
            params = target_string(params, params.global_idx, k_target_papillary, eta_papillary);     
        else
            params = target_string(params, params.global_idx, k_target_papillary);     
        end 
        
        % write papillary file 
        params = papillary_string(params, params.global_idx, leaflet.chordae(tree_idx).root); 
        
        params.global_idx = params.global_idx + 1;
        
        for i=1:N_chordae
            params.vertices(:,params.global_idx + 1) = C(:,i); 
            leaflet.chordae(tree_idx).indices_global(i) = params.global_idx;                 
            params.global_idx = params.global_idx + 1;
        end 
        
    end 

end 
 

function params = add_springs(params, leaflet, ds, collagen_spring)

    params = add_leaflet_springs(params, leaflet, ds, collagen_spring); 
    params = add_chordae_tree_springs(params, leaflet, ds, collagen_spring); 

end 


function params = add_leaflet_springs(params, leaflet, ds, collagen_spring)
                      
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    % Unpack needed data 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    is_bc             = leaflet.is_bc; 
    chordae           = leaflet.chordae;
    chordae_idx       = leaflet.chordae_idx; 

    R_u               = leaflet.R_u;
    k_u               = leaflet.k_u;
    R_v               = leaflet.R_v;
    k_v               = leaflet.k_v;

       
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
   
    % output flag information 
    copy = params.copy; 
    if params.output.leaflets(copy)
        output        = true; 
        output_stride = params.output.stride_leaflet; 
    else 
        output        = false; 
    end 
    
    if params.output.chordae(copy)
        output_tmp_chordae = true; 
    end 
    
    
    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point may have springs connected to it 
            if is_internal(j,k) || is_bc(j,k)
                
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp_k = true; 
                else 
                    output_tmp_k = false; 
                end 
                
                if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                    output_tmp_j = true; 
                else 
                    output_tmp_j = false; 
                end                 
                
                % global index of current point 
                idx = leaflet.indices_global(j,k); 
                
                % current node has a chordae connection
                if chordae_idx(j,k).tree_idx
                    
                    tree_idx = chordae_idx(j,k).tree_idx; 
                    
                    [m N_chordae] = size(chordae(tree_idx).C);

                    % index in current free edge array 
                    i = chordae_idx(j,k).leaf_idx;
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet
                    leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                    % then take the parent index of that number in chordae variables
                    idx_chordae = floor(leaf_idx/2);  
                    
                    nbr_idx = chordae(tree_idx).indices_global(idx_chordae); 

                    rest_len = chordae(tree_idx).R_free_edge(i); 

                    k_rel = chordae(tree_idx).k_free_edge(i); 
                    
                    params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output_tmp_chordae); 
                    
                end 
                                
                % springs in leaflet, only go in up direction 
                j_nbr_tmp = j + 1; 
                k_nbr_tmp = k; 
                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                if valid 
                    
                    % no bc to bc springs 
                    if ~(is_bc(j, k) && is_bc(j_nbr, k_nbr))
                    
                        rest_len = R_u(j_spr, k_spr); 
                        k_rel    = k_u(j_spr, k_spr); 

                        nbr_idx = leaflet.indices_global(j_nbr,k_nbr);
                        
                        
                        if j_nbr_tmp ~= j_nbr 
                            % periodic wrapping requires oppositite order 
                            params = place_spring_and_split(params, nbr_idx, idx, k_rel, rest_len, ds, collagen_spring, output_tmp_j);
                        else 
                           % standard order 
                            params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output_tmp_j);
                        end 

                    end 

                end 
                
                % springs in leaflet, only go in up direction 
                j_nbr_tmp = j; 
                k_nbr_tmp = k + 1; 
                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                if valid 
                    
                    % no bc to bc springs 
                    if ~(is_bc(j, k) && is_bc(j_nbr, k_nbr))
                    
                        rest_len = R_v(j_spr, k_spr); 
                        k_rel    = k_v(j_spr, k_spr); 

                        nbr_idx = leaflet.indices_global(j_nbr,k_nbr);

                        params = place_spring_and_split(params, idx, nbr_idx, k_rel, rest_len, ds, collagen_spring, output_tmp_k);

                    end 

                end 
                
            end 
        end 
    end
end 



function params = add_chordae_tree_springs(params, leaflet, ds, collagen_spring)
    % 
    % Adds chordae tree to IBAMR format files 
    % 
    
    num_trees = leaflet.num_trees; 
    chordae   = leaflet.chordae; 
        
    copy = params.copy; 
    if params.output.chordae(copy)
        output_tmp = true;
    else 
        output_tmp = false;
    end 
    
    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        [m N_chordae] = size(C);   
        
        indices_global = chordae(tree_idx).indices_global; 
        
        for i=1:N_chordae

            idx = indices_global(i); 
            
            % place the spring which goes to the parent 
            parent = floor(i/2); 
            
            if parent == 0
                nbr_idx = chordae(tree_idx).idx_root; 
            else
                nbr_idx = indices_global(parent); 
            end 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr rest_len k_rel] = get_nbr_chordae(leaflet, i, parent, tree_idx); 
            
            % list nbr index first because nbr is parent and has lower index
            params = place_spring_and_split(params, nbr_idx, idx, k_rel, rest_len, ds, collagen_spring, output_tmp);
                
        end 
        
    end 

end 


function params = place_cross_layer_springs(params)

    function_idx = 0; 
    kappa        = params.kappa_cross_layer; 
    rest_len     = params.rest_len_cross_layer;
    % don't include these for now 
    output_tmp = false; 
        
    for i=(params.min_idx_for_cross_layer):(params.max_idx_for_cross_layer)
        idx          = i -  params.total_per_layer; 
        nbr_idx      = i;  

        params = spring_string(params, idx, nbr_idx, kappa, rest_len, function_idx, output_tmp); 
    end 

end 


                        
function params = place_net(params, leaflet, ds, r, L, k_rel, k_target, ref_frac, eta)
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

    X     = leaflet.X; 
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
    
    % output flag information 
    copy = params.copy; 
    if params.output.mesh(copy)
        output        = true; 
        output_stride = params.output.stride_mesh; 
    else 
        output        = false; 
    end 
   
    function_idx = 0; 
    
    % only include those full rings which fit in the domain 
    % always place at least one, which is the ring itself 
    k_max_rings = max(1, floor( (L-r) / ds)); 
        
    points = zeros(3,j_max,k_max_rings);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(j_max,k_max_rings); 
    
    instrument_idx = 0; 
    
    % compute vertex positions and add to array 
    for k=1:k_max_rings
        for j=1:j_max
            
            % valve ring points from leaflet 
            ring_pt = X(1:2,j,k_max);

            increment = ring_pt - params.ring_center(1:2); 
            increment = ds * increment / norm(increment); 
            
            % expand in the direction of a vector from venter of ring to
            % current point 
            coords_horiz = ring_pt + (k-1)*increment; 
              
            % alternative, just scalar mutiply the point 
            % coords_horiz = (1 + (k-1)*ds) * ring_pt; 
                        
            % these might be the same... 
            
            % if one norm is less than L, then the point is within the domain  
            %if norm(coords_horiz, inf) < L 
            if (params.x_min    <= coords_horiz(1)) && ...   
               (coords_horiz(1) <= params.x_max   ) && ...   
               (params.y_min    <= coords_horiz(2)) && ...   
               (coords_horiz(2) <= params.y_max   ) 
                
                points(:,j,k) = [coords_horiz; params.ring_center(3)]; 
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

    
    % write the instrument file header here 
    fprintf(params.inst, '1   # num meters in file\n'); 
    fprintf(params.inst, 'meter_0   # name\n'); 
    fprintf(params.inst, '%d  # number of meter points\n', j_max); 
    
    % below the first possible point 
    idx = -1; 
    
    for k=1:k_max_rings
        for j=1:j_max
            
            % take mod one be 
            if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                output_tmp = true;
            else 
                output_tmp = false;
            end 
            
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
                if j < j_max
                    j_nbr = j+1; 
                    if ~isnan(indices_global(j_nbr,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j_nbr,k)); 
                        k_abs = k_rel / rest_len;
                        nbr_idx = indices_global(j_nbr,k); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp); 
                    end 
                end 
                
                % don't forget the periodic direction in j
                if j == j_max
                   j_nbr = 1; 
                   % need to make sure that the 1,k point is also not a NaN  
                   if ~isnan(indices_global(j_nbr,k)) 
                       rest_len = ref_frac * norm(points(:,j,k) - points(:,j_nbr,k)); 
                       k_abs = k_rel / rest_len;
                       nbr_idx = indices_global(j_nbr,k); 
                       params = spring_string(params, nbr_idx, idx, k_abs, rest_len, function_idx, output_tmp); 
                   end 
                end 
                
            end 
           
        end 
    end  

end 


function params = place_rays(params, leaflet, ds, L, k_rel, k_target, ref_frac, eta)
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

    X            = leaflet.X; 
    j_max        = leaflet.j_max; 
    k_max        = leaflet.k_max;
    is_bc        = leaflet.is_bc; 
    is_internal  = leaflet.is_internal; 
    
    function_idx = 0; 
    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
    
    if ~isfield(leaflet, 'indices_global')
        error('Must place leaflets before placing rays'); 
    end 
    
    % output flag information 
    copy = params.copy; 
    if params.output.mesh(copy)
        output        = true; 
        
        % use leaflet stride here so fibers that continue as rays are plotted as such
        output_stride = params.output.stride_leaflet; 
    else 
        output        = false; 
    end 
    
    
    for j = 1:j_max
        for k = 1:k_max
            if is_bc(j,k)
            
                % take mod one be 
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp = true;
                else 
                    output_tmp = false;
                end 
                
                pt_ring = X(:,j,k); 

                % only get a fiber if the previous point is included in the leaflet  
                neighbors = []; 
                
                % possible to have both directions of j nbr
                for j_nbr = [j-1,j+1]
                    k_nbr = k; 
                    if (j_nbr > 0) &&  (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && is_internal(j_nbr,k_nbr)
                        neighbors = [X(:,j_nbr,k_nbr), neighbors] ; 
                    end
                end 
                
                % k_nbr always down 
                j_nbr = j; 
                k_nbr = k-1; 
                if (j_nbr > 0) &&  (k_nbr > 0) && (j_nbr <= j_max) && (k_nbr <= k_max) && is_internal(j_nbr,k_nbr)
                    neighbors = [X(:,j_nbr,k_nbr), neighbors] ; 
                end        

                for x = neighbors 

                    % find the initial reflected point 
                    
                    % adjacent ring points determine local normal and tangent 
                    j_plus__1 = get_j_nbr(j+1, k, periodic_j, j_max); 
                    j_minus_1 = get_j_nbr(j-1, k, periodic_j, j_max);
                    
                    ring_nbr_plus  = X(:, j_plus__1, k); 
                    ring_nbr_minus = X(:, j_minus_1, k); 
                    
                    val = get_geodesic_continued_point(x, pt_ring, ring_nbr_minus, ring_nbr_plus, params.ring_center); 

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
                    % while norm(point(1:2), inf) < L   
                    while (params.x_min <= point(1))      && ...   
                       (point(1)     <= params.x_max)  && ...   
                       (params.y_min <= point(2))      && ...   
                       (point(2)     <= params.y_max) 
           
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

                        params = spring_string(params, nbr_idx, idx, k_abs, rest_len, function_idx, output_tmp);

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


function [val] = get_geodesic_continued_point(x, pt_ring, ring_nbr_minus, ring_nbr_plus, ring_center)
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
    
    if abs(pt_ring(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    if abs(ring_nbr_plus(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    if abs(ring_nbr_minus(3) - ring_center(3)) > tol
        error('Initial ring point is not near z ring plane'); 
    end 
    
    % local tangent implies local normal 
    tangent = ring_nbr_plus(1:2) - ring_nbr_minus(1:2); 
    normal  = [tangent(2); -tangent(1)]; 
    normal  = normal / norm(normal); 
    
    % translate ring point to origin (implicitly)
    % and x somewhere near the origin 
    val = x - pt_ring; 

    % rotate system around z such that normal now points in x axis direction 
    theta = atan2(normal(2), normal(1)); 
    
    val = rotation_matrix_z(-theta) * val; 
    
    % rotate val into the z = 0 plane, inside the transformed ring 
    % this is not inverted as it gets the image of the geodesic point before reflection 
    phi = atan2(val(3), val(1)); 
    val = rotation_matrix_y( -(phi - pi)) * val;
    
    if abs(val(3)) > tol
        error('should have landed in z=0 plane here...');
    end 

    % reflect, this would be the geodesic point if the system was flat 
    val = -val; 
    
    % undo rigid rotations and translations 
    val = rotation_matrix_z(theta) * val; 
    val = val + pt_ring; 
    
end 






function params = place_cartesian_net(params, leaflet, r_extra, L, ds, k_rel, k_target, ref_frac, eta)
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

    % compute the valve ring, adding r_extra to each point 
    X     = leaflet.X; 
    j_max = leaflet.j_max; 
    k_max = leaflet.k_max; 
   
    function_idx = 0; 
    
    ring_expanded = zeros(2,j_max); 
    
    for j=1:j_max
        % valve ring points from leaflet 
%         ring_pt    = X(1:2,j,k_max); 
%         increment  = r_extra * ring_pt / norm(ring_pt); 
% 
%         ring_expanded(:,j) = ring_pt + increment; 

        % valve ring points from leaflet 
        ring_pt = X(1:2,j,k_max);

        increment = ring_pt - params.ring_center(1:2); 
        increment = r_extra * increment / norm(increment); 

        % expand in the direction of a vector from venter of ring to
        % current point 
        ring_expanded(:,j) = ring_pt + increment; 

    end 
    
    % output flag information 
    copy = params.copy; 
    if params.output.cartesian_mesh(copy)
        output        = true; 
        output_stride = params.output.stride_mesh; 
    else 
        output        = false; 
    end 
    
    
    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices_global = zeros(N,N);     
    
    
    % This loop should be one indexed, 
    % because we want to start not at the edge but ds in 
    for k=1:N
        for j=1:N
                        
            % coords_horiz = [ (j-1)*ds - L + ds/2; (k-1)*ds - L + ds/2]; 
            coords_horiz = [ (j-1)*ds + params.x_min + ds/2; (k-1)*ds + params.y_min + ds/2]; 
            
            % if one norm is less than L, then the point is within the domain  
            % if (norm(coords_horiz, inf) < L) && point_out_of_polygon(ring_expanded, coords_horiz)
            if (params.x_min    <= coords_horiz(1)) && ...   
               (coords_horiz(1) <= params.x_max   ) && ...   
               (params.y_min    <= coords_horiz(2)) && ...   
               (coords_horiz(2) <= params.y_max   ) && ... 
               point_out_of_polygon(ring_expanded, coords_horiz)
                
                points(:,j,k) = [coords_horiz; params.ring_center(3)]; 
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
 
                if output && ((mod(j,output_stride) == 1) || (output_stride == 1))
                    output_tmp_k = true; 
                else 
                    output_tmp_k = false; 
                end 
                
                if output && ((mod(k,output_stride) == 1) || (output_stride == 1))
                    output_tmp_j = true; 
                else 
                    output_tmp_j = false; 
                end    
                
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
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp_j); 
                    end 
                end 
                                
                % check up directions for springs 
                if (k+1) <= N
                    if ~isnan(indices_global(j,k+1))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                        k_abs = k_rel / rest_len; 
                        nbr_idx = indices_global(j,k+1); 
                        params = spring_string(params, idx, nbr_idx, k_abs, rest_len, function_idx, output_tmp_k); 
                    end 
                end 
                
            end 
           
        end 
    end 

end 


function outside = point_out_of_polygon(vert, test)
%
% Recklessly translated from 
% 
% https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
%
% Original code: 
% 
% int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
% {
%   int i, j, c = 0;
%   for (i = 0, j = nvert-1; i < nvert; j = i++) {
%     if ( ((verty[i]>testy) != (verty[j]>testy)) &&
% 	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
%        c = !c;
%   }
%   return c;
% }

    inside = false; 
    
    if size(vert,1) ~= 2
        error('Must pass two d vector of vertices to test');         
    end 
    
    nvert = size(vert,2); 
    
    i = 1; 
    j = nvert; 
    
    while i <= nvert 
        
        check_1 = ((vert(2,i) > test(2)) ~= (vert(2,j) > test(2))); 
        tmp     = (vert(1,j) - vert(1,i)) * (test(2)-vert(2,i)) / (vert(2,j) - vert(2,i)) + vert(1,i); 
        check_2 = (test(1) < tmp); 
        
        if check_1 && check_2
            inside = ~inside; 
        end 
        
        j = i; 
        i = i+1; 
    end 
    
    outside = ~inside; 

end 


function [params, total_lagrangian_placed] = place_lagrangian_tracers(params, n_lagrangian_tracers, double_z)
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
    
    % includes ring_center
    % do not need to update this for nonzero ring_center
    x_min = params.x_min; 
    x_max = params.x_max; 
    y_min = params.y_min; 
    y_max = params.y_max; 
    z_min = params.z_min; 
    z_max = params.z_max;
    
    n_z_dir = n_lagrangian_tracers;  
    if double_z
        n_z_dir = 2*n_z_dir; 
    end 
    
    dx = (x_max - x_min) / n_lagrangian_tracers; 
    dy = (y_max - y_min) / n_lagrangian_tracers; 
    dz = (z_max - z_min) / n_z_dir; 
    
    total_lagrangian_placed = 0; 
    
    for i = 1:n_lagrangian_tracers
        for j = 1:n_lagrangian_tracers
            for k = 1:n_z_dir 
                
                x = (i - .5)*dx + x_min; 
                y = (j - .5)*dy + y_min; 
                z = (k - .5)*dz + z_min; 
                
                params.vertices(:,params.global_idx + 1) = [x y z]; 
    
                total_lagrangian_placed = total_lagrangian_placed + 1; 
                params.global_idx = params.global_idx + 1; 
            end 
        end 
    end 

end 





