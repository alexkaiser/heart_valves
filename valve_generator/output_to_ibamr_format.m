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
    
    N                            = valve.N; 
    base_name                    = valve.base_name; 
    L                            = valve.L; 
    pressure_tension_ratio       = valve.pressure_tension_ratio; 
    posterior                    = valve.posterior; 
    anterior                     = valve.anterior; 
    p_physical                   = valve.p_physical; 
    target_multiplier            = valve.target_multiplier; 
    refinement                   = valve.refinement; 
    n_lagrangian_tracers         = valve.n_lagrangian_tracers; 
    num_copies                   = valve.num_copies; 
    collagen_springs_leaflet     = valve.collagen_springs_leaflet; 
    

    vertex = fopen(strcat(base_name, '.vertex'), 'w'); 
    spring = fopen(strcat(base_name, '.spring'), 'w'); 
    target = fopen(strcat(base_name, '.target'), 'w'); 
    inst   = fopen(strcat(base_name, '.inst'  ), 'w'); 

    % keep one global index through the whole thing 
    global_idx = 0;
    
    % just count the number of vertices and strings throughout 
    total_vertices = 0; 
    total_springs  = 0; 
    total_targets  = 0; 

    % compute some needed constants 
    MMHG_TO_CGS = 1333.22368; 
    p_cgs = p_physical * MMHG_TO_CGS; 

    % value of the relative spring constant is determined by the ratio 
    k_rel = p_cgs / pressure_tension_ratio; 
    k_rel = k_rel / num_copies; 
    
    % base rate for target spring constants
    % target constant for a single point 
    % this does not scale when the mesh is changed 
    k_target = target_multiplier * k_rel; 
    
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
    
    % relative spring constants drop when the mesh is refined 
    k_rel = k_rel / refinement; 

    % output the left and right papillary as the first two vertices and targets
    
    % Critical damping for given k, mass m is 2*sqrt(m*k) 
    % Set to half critical for first test
    % 
    % m_effective_papillary = 1.0 * pi * filter_params_posterior.r^2; 
    eta_papillary         = 0.0; %sqrt(k_target/2 * m_effective_papillary); 
    
    % spacing 
    ds = 2*L / N;
    
    k_rel_leaflet = k_rel; 
    if collagen_springs_leaflet
        k_rel_leaflet = ds / num_copies; 
    end 
    
    if posterior.reflect_x
        posterior.X(1,:,:)                   = -posterior.X(1,:,:); 
        posterior.R(1,:,:)                   = -posterior.R(1,:,:); 
        posterior.left_papillary(1)          = -posterior.left_papillary(1); 
        posterior.right_papillary(1)         = -posterior.right_papillary(1); 
        posterior.chordae.C_left (1,:,:)     = -posterior.chordae.C_left (1,:,:);
        posterior.chordae.C_right(1,:,:)     = -posterior.chordae.C_right(1,:,:);
        posterior.chordae.left_papillary(1)  = -posterior.chordae.left_papillary(1); 
        posterior.chordae.right_papillary(1) = -posterior.chordae.right_papillary(1); 
    end 
    
    if anterior.reflect_x
        warning('Something strange, should not be reflecting on anterior leaflet')
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
    
    
    if valve.X_config_is_reference
        anterior.R               = anterior.X; 
        anterior.chordae.Ref_l   = anterior.chordae.C_left;  
        anterior.chordae.Ref_r   = anterior.chordae.C_right;  
        anterior.ref_frac        = 1.0; 
        
        posterior.R              = posterior.X; 
        posterior.chordae.Ref_l  = posterior.chordae.C_left;  
        posterior.chordae.Ref_r  = posterior.chordae.C_right;  
        posterior.ref_frac       = 1.0;
    end 
    
    % ugh this is terrible fix it it makes my head hurt by I'm tired 
    global z_offset
    
    for z_offset = z_offset_vals
        
        left_papillary  = posterior.left_papillary; 
        posterior.left_papillary_idx = global_idx; 
        total_vertices  = vertex_string(vertex, left_papillary, total_vertices); 
        total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
        global_idx      = global_idx + 1; 

        right_papillary = posterior.right_papillary; 
        posterior.right_papillary_idx = global_idx; 
        total_vertices  = vertex_string(vertex, right_papillary, total_vertices); 
        total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
        global_idx      = global_idx + 1;     

        if valve.split_papillary

            left_papillary  = anterior.left_papillary; 
            anterior.left_papillary_idx = global_idx; 
            total_vertices  = vertex_string(vertex, left_papillary, total_vertices); 
            total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
            global_idx      = global_idx + 1; 

            right_papillary = anterior.right_papillary; 
            anterior.right_papillary_idx = global_idx; 
            total_vertices  = vertex_string(vertex, right_papillary, total_vertices); 
            total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
            global_idx      = global_idx + 1;             

        else
            anterior.left_papillary_idx  = posterior.left_papillary_idx;
            anterior.right_papillary_idx = posterior.right_papillary_idx;  
        end 

        
        % posterior first 
        % leaflets 
        [global_idx, total_vertices, total_springs, total_targets, posterior] = ...
            add_leaflet(posterior, spring, vertex, target, ...
                        global_idx, total_vertices, total_springs, total_targets, k_rel_leaflet, k_target_ring, eta_ring, collagen_springs_leaflet); 

        % posterior chordae 
        posterior.chordae.left_papillary_idx  = posterior.left_papillary_idx; 
        posterior.chordae.right_papillary_idx = posterior.right_papillary_idx; 

        [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(posterior, spring, vertex, global_idx, total_vertices, total_springs, k_rel_leaflet, collagen_springs_leaflet);  

        % anterior 
        [global_idx, total_vertices, total_springs, total_targets, anterior] = ...
            add_leaflet(anterior, spring, vertex, target, ...
                         global_idx, total_vertices, total_springs, total_targets, k_rel_leaflet, k_target_ring, eta_ring, collagen_springs_leaflet);   


        % anterior chordae 
        anterior.chordae.left_papillary_idx  = anterior.left_papillary_idx; 
        anterior.chordae.right_papillary_idx = anterior.right_papillary_idx; 

        [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(anterior, spring, vertex, global_idx, total_vertices, total_springs, k_rel_leaflet, collagen_springs_leaflet);  

        % flat part of mesh 
        r = valve.r; 
        N_ring = 2 * N;
        h = 0.0; % ring always set at zero 
        ref_frac_net = 1.0; 

        % no radial fibers, instead geodesics from the leaflet 
        radial_fibers = false; 

        % turn the polar net off for now      
        [global_idx, total_vertices, total_springs, total_targets] = ...
                                place_net(r, h, L, N_ring, radial_fibers, spring, vertex, target, inst, ...
                                global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net); 

                            
        % place rays for now 
        [global_idx, total_vertices, total_springs, total_targets] = ...
                                place_rays(anterior, L, spring, vertex, target, ...
                                        global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net);                    

        [global_idx, total_vertices, total_springs, total_targets] = ...
                                place_rays(posterior, L, spring, vertex, target, ...
                                        global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net);                    


        % flat part of mesh with Cartesian coordinates
        % inner radius, stop mesh here 
        r_cartesian = r + 2*ds; 
        [global_idx, total_vertices, total_springs, total_targets] = ...
                                place_cartesian_net(r_cartesian, h, L, ds, spring, vertex, target, ...
                                global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net); 
 
    end 
                        
                        
    if n_lagrangian_tracers > 0
        double_z = false; 
        [global_idx, total_vertices, total_lagrangian_placed] = place_lagrangian_tracers(global_idx, total_vertices, vertex, n_lagrangian_tracers, L, double_z); 
        particles = fopen(strcat(base_name, '.particles'), 'w'); 
        fprintf(particles, '%d\n' ,total_lagrangian_placed); 
    end 

    % clean up files with totals 
    fclose(vertex); 
    fclose(spring); 
    fclose(target); 
    fclose(inst  ); 

    prepend_line_with_int(strcat(base_name, '.vertex'), total_vertices); 
    prepend_line_with_int(strcat(base_name, '.spring'), total_springs); 
    prepend_line_with_int(strcat(base_name, '.target'), total_targets); 

end 



% nest this function so it can access the z increment
% this is bad practice and should be removed 
function total_vertices = vertex_string(vertex, coords, total_vertices)
    % prints formatted string for current vertex to vertex file   
    
    % FIXME!!! 
    global z_offset
    
    fprintf(vertex, '%.14f\t %.14f\t %.14f\n', coords(1), coords(2), coords(3) + z_offset); 
    total_vertices = total_vertices + 1; 
end


function total_springs = spring_string(spring, idx, nbr, kappa, rest_len, total_springs, function_idx)
    % prints a spring format string to string file 
    if nbr <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    fprintf(spring, '%d\t %d\t %.14f\t %.14f', idx, nbr, kappa, rest_len); 
    
    % index for custom spring functions 
    if exist('function_idx', 'var') 
        fprintf(spring, '\t%d', function_idx); 
    end 
    
    fprintf(spring, '\n'); 
    
    total_springs = total_springs + 1; 
end 

function total_targets = target_string(target, idx, kappa, total_targets, eta)
    % prints a target format string to target file 
    
    if exist('eta', 'var') && (eta > 0.0)
        fprintf(target, '%d\t %.14f\t %.14f\n', idx, kappa, eta);
    else
        fprintf(target, '%d\t %.14f\n', idx, kappa);
    end 
    total_targets = total_targets + 1; 
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


function [global_idx, total_vertices, total_springs, total_targets, leaflet] = ...
                add_leaflet(leaflet, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, eta, collagen_spring)
                      
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    % Unpack needed data 
    X                 = leaflet.X; 
    R                 = leaflet.R; 
    is_internal       = leaflet.is_internal;
    point_idx_with_bc = leaflet.point_idx_with_bc; 
    is_bc             = leaflet.is_bc; 
    chordae           = leaflet.chordae;
    chordae_idx_left  = leaflet.chordae_idx_left; 
    chordae_idx_right = leaflet.chordae_idx_right; 
    ref_frac          = leaflet.ref_frac; 
    
    if collagen_spring
        function_idx = 1;
    end 
    
    % On some layouts N may not be the dimension 
    % because of boundary conditions 
    j_max = size(X,2); 
    k_max = size(X,2); 
    
    % Keep track of indices     
    indices_global = nan * zeros(j_max, k_max); 

    % Counter for number of points put down 
    pts_placed = 0; 
    
    if isfield(leaflet, 'chordae') && ~isempty(chordae)
        [m N_chordae] = size(leaflet.chordae.C_left); 
        total_leaflet = sum(is_internal(:)) + sum(is_bc(:)); 
    else 
        error('Leaflet without chordae not implemented'); 
    end 
    
    
    for k=1:k_max
        for j=1:j_max
            
            % every internal and boundary point written to the file 
            if is_internal(j,k) || is_bc(j,k)
                
                idx = global_idx + point_idx_with_bc(j,k);    
                
                total_vertices = vertex_string(vertex, X(:,j,k), total_vertices); 
                
                indices_global(j,k) = idx; 
                
                % Connect to the left papillary or chordae tree 
                if chordae_idx_left(j,k)
                    
                    % current node has a chordae connection
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_left(j,k) + N_chordae; 

                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  

                    R_nbr = leaflet.chordae.Ref_l(:,idx_chordae); 

                    nbr_idx = global_idx + total_leaflet + idx_chordae - 1; 

                    rest_len = ref_frac * norm(R_nbr - R(:,j,k)); 

                    if collagen_spring 
                        % relative constant here
                        % other parameters coded into function 
                        kappa = k_rel * chordae.k_0; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs, function_idx); 
                    else 
                        kappa = k_rel * chordae.k_0 / rest_len; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 

                end 
                
                % Connect to the right papillary or chordae tree 
                if chordae_idx_right(j,k)
                    
                    % current node has a chordae connection
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet 
                    leaf_idx = chordae_idx_right(j,k) + N_chordae; 

                    % then take the parent index of that number in chordae variables 
                    idx_chordae = floor(leaf_idx/2);  

                    R_nbr = leaflet.chordae.Ref_r(:,idx_chordae); 

                    nbr_idx = global_idx + total_leaflet + idx_chordae - 1; 
                    
                    % right gets incremented by N_chordae again     
                    nbr_idx = nbr_idx + N_chordae;     
                        
                    rest_len = ref_frac * norm(R_nbr - R(:,j,k)); 

                    if collagen_spring 
                        % relative constant here
                        % other parameters coded into function 
                        kappa = k_rel * chordae.k_0; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs, function_idx); 
                    else 
                        kappa = k_rel * chordae.k_0 / rest_len; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 

                end 
                               
                % if on boundary, this is a target point 
                if is_bc(j,k)
                    if exist('eta', 'var')
                        total_targets = target_string(target, idx, k_target, total_targets, eta);     
                    else
                        total_targets = target_string(target, idx, k_target, total_targets);     
                    end 
                end 
                
                % springs in leaflet, only go in up direction 
                j_nbr = j + 1; 
                k_nbr = k; 
                if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))
                    
                    rest_len = ref_frac * norm(R(:,j_nbr,k_nbr) - R(:,j,k)); 
                    nbr_idx = global_idx + point_idx_with_bc(j_nbr,k_nbr);
                    
                    if collagen_spring
                        kappa = k_rel;         
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs, function_idx); 
                    else 
                        kappa = k_rel / rest_len;         
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 

                end 
                
                % springs in leaflet, only go in up direction 
                j_nbr = j; 
                k_nbr = k + 1; 
                if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))
                    
                    rest_len = ref_frac * norm(R(:,j_nbr,k_nbr) - R(:,j,k)); 
                    nbr_idx = global_idx + point_idx_with_bc(j_nbr,k_nbr);
                    
                    if collagen_spring
                        kappa = k_rel;         
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs, function_idx); 
                    else 
                        kappa = k_rel / rest_len;         
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 

                end 
                
                pts_placed = pts_placed + 1; 
            end 
        end 
    end
    
    leaflet.indices_global = indices_global; 
    global_idx = global_idx + pts_placed;

end 



function [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(leaflet, spring, vertex, global_idx, total_vertices, total_springs, k_rel, collagen_spring)

    % Adds chordae tree to IBAMR format files 
    % No targets here, so files and and count not included 
                        
    if ~isfield(leaflet, 'chordae') || isempty(leaflet.chordae)
        error('cannot place chordae with empty array'); 
    end 
    
    chordae         = leaflet.chordae; 
    ref_frac        = leaflet.ref_frac; 
    
    C_left          = chordae.C_left;  
    C_right         = chordae.C_right;  
    Ref_l           = chordae.Ref_l;  
    Ref_r           = chordae.Ref_r; 

    
    [m N_chordae] = size(chordae.C_left); 
    
    if collagen_spring
        function_idx = 1;
    end
    
    
    for left_side = [true false];  
        
        if left_side
            C = C_left; 
            Ref = Ref_l; 
        else 
            C = C_right; 
            Ref = Ref_r; 
        end 
        
        for i=1:N_chordae

            idx = global_idx + i - 1; % subtract one ibamr is zero indexed 
            
            % right side gets an extra factor of N_chordae
            if ~left_side
                idx = idx + N_chordae; 
            end 
            
            % the vertex is always placed 
            total_vertices = vertex_string(vertex, C(:,i), total_vertices); 
            
            % place the spring which goes to the parent 
            parent = floor(i/2); 
            
            if parent == 0
                
                % zero index means papillary muscles
                if left_side 
                    nbr_idx = chordae.left_papillary_idx; 
                else 
                    nbr_idx = chordae.right_papillary_idx; 
                end
                
            else
                nbr_idx = parent + global_idx - 1; % subtract one, ibamr zero indexed  
                if ~left_side
                    nbr_idx = nbr_idx + N_chordae; 
                end 
            end 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, parent, left_side); 
            
            rest_len = ref_frac * norm(R_nbr - Ref(:,i)); 
            
            if collagen_spring 
                kappa = k_rel * k_val; 
                % list nbr index first because nbr is parent and has lower index
                total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs, function_idx);            
            else 
                kappa = k_rel * k_val / rest_len; 
                % list nbr index first because nbr is parent and has lower index
                total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs);        
            end 
        end 
        
    end 
    
    global_idx = global_idx + 2*N_chordae; 

end 
                        
                        
function [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_net(r, h, L, N, radial_fibers, spring, vertex, target, inst, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, ref_frac, eta)
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
    %     global_idx       Running totals  
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
    indices = zeros(N,M); 
    idx = 0; 
    points_placed = 0; 
    

    rad = r; 
    
    for k=1:M
        for j=1:N
            
            theta = (j-1) * ds; 
            coords_horiz = [rad*cos(theta); rad*sin(theta)]; 
            
            % if one norm is less than L, then the point is within the domain  
            if norm(coords_horiz, inf) < L 
                points(:,j,k) = [coords_horiz; h]; 
                indices(j,k) = idx; 
                idx = idx + 1;                 
            else 
                points(:,j,k) = NaN * ones(3,1);
                indices(j,k) = NaN; 
            end     
            
        end 
        
        rad = rad + ds; 
    end 

    
    % write the instrument file header here 
    fprintf(inst, '1   # num meters in file\n'); 
    fprintf(inst, 'meter_0   # name\n'); 
    fprintf(inst, '%d  # number of meter points\n', N); 
    
    
    % below the first possible point 
    idx = -1; 
    
    for k=1:M
        for j=1:N
            
            % just ignore the nan 
            if ~isnan(indices(j,k))
 
                last_idx = idx; 
                idx = indices(j,k) + global_idx; 
                
                % instrument file on 
                if k == 1
                
                    if j ~= (indices(j,k)+1)
                        error('local indices should go with j'); 
                    end 
                
                    fprintf(inst, '%d \t0 \t %d\n', idx, indices(j,k)); 
                    
                end 
                    
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 
                
                % every valid vertex is a target here 
                total_vertices = vertex_string(vertex, points(:,j,k), total_vertices); 
                points_placed = points_placed + 1; 
                if exist('eta', 'var')
                    total_targets = target_string(target, idx, k_target, total_targets, eta);     
                else
                    total_targets = target_string(target, idx, k_target, total_targets);     
                end 
                
                
                % check up directions for springs 
                if (j+1) < N
                    if ~isnan(indices(j+1,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j+1,k)); 
                        kappa = k_rel / rest_len;
                        nbr_idx = indices(j+1,k) + global_idx; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 
                end 
                
                % don't forget the periodic direction in j
                if (j+1) == N
                   % need to make sure that the 1,k point is also not a NaN  
                   if ~isnan(indices(1,k)) 
                       rest_len = ref_frac * norm(points(:,j,k) - points(:,1,k)); 
                       kappa = k_rel / rest_len;
                       nbr_idx = indices(1,k) + global_idx; 
                       total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs); 
                   end 
                end 
                
                % check up directions for springs 
                if radial_fibers
                    if (k+1) < M
                        if ~isnan(indices(j,k+1))
                            rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                            kappa = k_rel / rest_len; 
                            nbr_idx = indices(j,k+1) + global_idx; 
                            total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                        end 
                    end 
                end 
                
                % no periodic direction in radial direction (k)
                
            end 
           
        end 
    end 
    
    global_idx = global_idx + points_placed; 

end 


function [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_rays(leaflet, L, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, ref_frac, eta)
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
    %     global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    X           = leaflet.X; 
    is_bc       = leaflet.is_bc; 
    is_internal = leaflet.is_internal; 
    r           = leaflet.filter.r; 
    h           = 0.0;        % always place at origin 
    N           = leaflet.N; 
    
    if ~isfield(leaflet, 'indices_global')
        error('Must place leaflets before placing rays'); 
    end 
    
    j_max = size(X,2); 
    k_max = size(X,2); 
    
    for j = 1:j_max
        for k=1:k_max
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


                    point = val; 
                    point_prev = pt_ring; 
                    nbr_idx = leaflet.indices_global(j,k); 

                    % just keep adding until points leave the domain 
                    while norm(point(1:2), inf) < L   

                        % grab the index 
                        idx = global_idx;

                        % point 
                        total_vertices = vertex_string(vertex, point, total_vertices); 

                        % it's a target too 
                        if exist('eta', 'var')
                            total_targets = target_string(target, idx, k_target, total_targets, eta);     
                        else
                            total_targets = target_string(target, idx, k_target, total_targets);     
                        end 


                        rest_len = ref_frac * norm(point - point_prev); 
                        kappa = k_rel / rest_len;


                        total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs);

                        point_prev = point; 
                        point      = point + increment; 
                        global_idx = global_idx + 1; 
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






function [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_cartesian_net(r, h, L, ds, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, ref_frac, eta)
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
    %     global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    
    N = ceil(2*L / ds);  
    points = zeros(3,N,N);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices = zeros(N,N); 
    idx = 0; 
    points_placed = 0; 
    
    
    % This loop should be one indexed, 
    % because we want to start not at the edge but ds in 
    for k=1:N
        for j=1:N
            
            coords_horiz = [ (j-1)*ds - L + ds/2; (k-1)*ds - L + ds/2]; 
            
            % if one norm is less than L, then the point is within the domain  
            if (norm(coords_horiz, inf) < L) && (norm(coords_horiz,2) > r) 
                points(:,j,k) = [coords_horiz; h]; 
                indices(j,k) = idx; 
                idx = idx + 1;                 
            else 
                points(:,j,k) = NaN * ones(3,1);
                indices(j,k) = NaN; 
            end     
            
        end 
        
    end 

    
    % below the first possible point 
    idx = -1; 
    
    for k=1:N
        for j=1:N
            
            % just ignore the NaNs 
            if ~isnan(indices(j,k))
 
                last_idx = idx; 
                idx = indices(j,k) + global_idx; 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 
                
                % every valid vertex is a target here 
                total_vertices = vertex_string(vertex, points(:,j,k), total_vertices); 
                points_placed = points_placed + 1; 
                if exist('eta', 'var')
                    total_targets = target_string(target, idx, k_target, total_targets, eta);     
                else
                    total_targets = target_string(target, idx, k_target, total_targets);     
                end 

                
                % check up directions for springs 
                if (j+1) <= N
                    if ~isnan(indices(j+1,k))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j+1,k)); 
                        kappa = k_rel / rest_len;
                        nbr_idx = indices(j+1,k) + global_idx; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 
                end 
                                
                % check up directions for springs 
                if (k+1) <= N
                    if ~isnan(indices(j,k+1))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                        kappa = k_rel / rest_len; 
                        nbr_idx = indices(j,k+1) + global_idx; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 
                end 
                
            end 
           
        end 
    end 
    
    global_idx = global_idx + points_placed; 

end 



function [global_idx, total_vertices, total_lagrangian_placed] = place_lagrangian_tracers(global_idx, total_vertices, vertex, n_lagrangian_tracers, L, double_z)
    % Places a uniform cartesian mesh of lagrangian particle tracers 
    % Simple lopp implementation 
    %
    %     global_idx               Running totals  
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
                
                total_vertices = vertex_string(vertex, [x y z], total_vertices); 
    
                total_lagrangian_placed = total_lagrangian_placed + 1; 
                global_idx = global_idx + 1; 
            end 
        end 
    end 

end 





