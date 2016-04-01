function [] = output_to_ibamr_format(base_name, L, ratio, params_posterior, filter_params_posterior, params_anterior, p_physical, target_multiplier, refinement, n_lagrangian_tracers)
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
    % 
    % Output: 
    %    Files written in IBAMR format 
    % 

    vertex = fopen(strcat(base_name, '.vertex'), 'w'); 
    spring = fopen(strcat(base_name, '.spring'), 'w'); 
    target = fopen(strcat(base_name, '.target'), 'w'); 
    inst = fopen(strcat(base_name, '.inst'  ), 'w'); 

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
    k_rel = p_cgs / ratio; 
    
    % base rate for target spring constants
    % target constant for a single point 
    % this does not scale when the mesh is changed 
    k_target = target_multiplier * k_rel; 
    
    % the valve ring is 1d, should be halfed with doubling of mesh 
    % also set damping coefficients accordingly 
    k_target_ring = k_target; %  / refinement; 
    m_ring = 0.001; 
    eta_ring = sqrt(m_ring * k_target_ring); 
    
    % there are four times as many, so they get multiplied by refinement squared 
    % can also just divide by refinement because not want them to get stiffer
    k_target_net = k_target; % / refinement; 
    m_net = 0.001; 
    eta_net = sqrt(m_net * k_target_net);    
    
    % relative spring constants drop when the mesh is refined 
    k_rel = k_rel / refinement; 

    % output the left and right papillary as the first two vertices and targets
    
    % Critical damping for given k, mass m is 2*sqrt(m*k) 
    % Set to half critical for first test 
    m_effective_papillary = .3; 
    eta_papillary         = sqrt(k_target * m_effective_papillary); 
    
    left_papillary  = [0; -filter_params_posterior.a; 0]; 
    total_vertices  = vertex_string(vertex, left_papillary, total_vertices); 
    total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
    global_idx      = global_idx + 1; 
    
    right_papillary = [0;  filter_params_posterior.a; 0]; 
    total_vertices  = vertex_string(vertex, right_papillary, total_vertices); 
    total_targets   = target_string(target, global_idx, k_target, total_targets, eta_papillary);     
    global_idx      = global_idx + 1;     
    
    % posterior first 
    % leaflets 
    [global_idx, total_vertices, total_springs, total_targets, params_posterior] = ...
        add_leaflet(params_posterior, filter_params_posterior, spring, vertex, target, ...
                    global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_ring, eta_ring); 

    % if chordae exist, then add them 
    if isfield(params_posterior, 'chordae') && ~isempty(params_posterior.chordae)
        [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(params_posterior, spring, vertex, global_idx, total_vertices, total_springs, k_rel);  
    end
    
    % anterior 
    [global_idx, total_vertices, total_springs, total_targets, params_anterior] = ...
        add_leaflet(params_anterior, filter_params_posterior, spring, vertex, target, ...
                     global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_ring, eta_ring);   
    
    if isfield(params_anterior, 'chordae') && ~isempty(params_anterior.chordae)
        [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(params_anterior, spring, vertex, global_idx, total_vertices, total_springs, k_rel);  
    end 
    
    % flat part of mesh 
    r = filter_params_posterior.r; 
    N_ring = 2 * params_anterior.N;
    h = filter_params_posterior.h; 
    ref_frac_net = 1.0; 
    
    % no radial fibers, instead geodesics from the leaflet 
    radial_fibers = false; 
    
    % turn the polar net off for now      
    [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_net(r, h, L, N_ring, radial_fibers, spring, vertex, target, inst, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net); 
                        
                        
    
    % place rays for now 
    [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_rays(params_anterior, filter_params_posterior, L, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net);                    
    
    [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_rays(params_posterior, filter_params_posterior, L, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net);                    
        
                        
    % flat part of mesh with Cartesian coordinates
    ds = 2*L / params_anterior.N;
    
    % inner radius, stop mesh here 
    r_cartesian = r + 2*ds; 
    [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_cartesian_net(r_cartesian, h, L, ds, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target_net, ref_frac_net, eta_net); 
 
                        
    if nargin >= 10
        double_z = false; 
        [global_idx, total_vertices, total_lagrangian_placed] = place_lagrangian_tracers(global_idx, total_vertices, vertex, n_lagrangian_tracers, L, filter_params_posterior, double_z); 
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



function total_vertices = vertex_string(vertex, coords, total_vertices)
    % prints formatted string for current vertex to vertex file   
    fprintf(vertex, '%.14f\t %.14f\t %14f\n', coords(1), coords(2), coords(3)); 
    total_vertices = total_vertices + 1; 
end 


function total_springs = spring_string(spring, idx, nbr, kappa, rest_len, total_springs)
    % prints a spring format string to string file 
    if nbr <= idx
        error('By convention, only place springs with the second index larger to prevent duplicates'); 
    end 
    
    fprintf(spring, '%d\t %d\t %.14f\t %.14f\n', idx, nbr, kappa, rest_len); 
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
    % 
    % Out of place and wasteful

    read = fopen(file_name, 'r'); 
    
    write = fopen(strcat(file_name, '.tmp'), 'w'); 
    
    % write the line...
    fprintf(write, '%d\n', val); 
    
    % now just copy all the lines over 
    tline = fgetl(read);
    while ischar(tline)
        fprintf(write, strcat(tline, '\n')); 
        tline = fgetl(read);
    end

    fclose(read); 
    fclose(write); 
    
    % clobber the old file with the temp file
    movefile(strcat(file_name, '.tmp'), file_name); 

end 

function [global_idx, total_vertices, total_springs, total_targets, params] = ...
                add_leaflet(params, filter_params, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, eta)
                                                                                
% Places all main data into IBAMR format for this leaflet
% Updates running totals on the way 


    [X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params); 
    
    % Keep track of indices 
    indices_global = nan * zeros(N+1,N+1); 

    
    left_papillary  = [0; -filter_params.a; 0]; 
    right_papillary = [0;  filter_params.a; 0]; 
    
    pts_placed = 0; 
    

    if isfield(params, 'chordae') && ~isempty(params.chordae)
        chordae_tree = true;
        [m N_chordae] = size(params.chordae.C_left); 
        total_leaflet = (N+1)*(N+2)/2; 
    else 
        chordae_tree = false;
    end 
    
    
    for k=1:N+1
        for j=1:N+1
            if (j+k) <= (N+2)
                
                idx = global_idx + vertex_index_offset(j,k,N);    
                
                total_vertices = vertex_string(vertex, X(:,j,k), total_vertices); 
                
                indices_global(j,k) = idx; 
                
                % if j==1, connect to the left papillary or chordae tree 
                if j==1
                    
                    if chordae_tree
                        
                        % there is one point which is a b.c. which is not included 
                        if is_internal(j,k,N)
                        
                            j_nbr = 0; 
                            k_nbr = k; 
                            [X_nbr R_nbr idx_chordae left_side] = get_neighbor(params, filter_params, j_nbr, k_nbr); 

                            nbr_idx = global_idx + total_leaflet + idx_chordae - 1; 

                            % right gets incremented by N_chordae again 
                            if ~left_side
                                nbr_idx = nbr_idx + N_chordae; 
                            end 

                            rest_len = ref_frac * norm(R_nbr - R(:,j,k)); 
                            kappa = k_rel * params.chordae.k_0 / rest_len; 

                            total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                        
                        end 
                    
                    else 
                        nbr_idx = 0; 
                        rest_len = ref_frac * norm(left_papillary - R(:,j,k)); 
                        kappa = k_rel / rest_len; 

                        total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs); 
                    end 
                end 
                
                % if k==1, connect to the right papillary or chordae tree 
                if k==1
                    
                    if chordae_tree
                        
                        % there is one point which is a b.c. which is not included 
                        if is_internal(j,k,N)
                        
                            j_nbr = j; 
                            k_nbr = 0; 
                            [X_nbr R_nbr idx_chordae left_side] = get_neighbor(params, filter_params, j_nbr, k_nbr); 

                            nbr_idx = global_idx + total_leaflet + idx_chordae - 1; 

                            % right gets incremented by N_chordae again 
                            if ~left_side
                                nbr_idx = nbr_idx + N_chordae; 
                            end 

                            rest_len = ref_frac * norm(R_nbr - R(:,j,k)); 
                            kappa = k_rel * params.chordae.k_0 / rest_len; 

                            total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                        
                        end 
                    else
                    
                        nbr_idx = 1; 
                        rest_len = ref_frac * norm(right_papillary - R(:,j,k)); 
                        kappa = k_rel / rest_len; 
                    
                        total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs); 
                    end 
                end                 
                
                % if on boundary, this is a target point 
                if (j+k) == (N+2)
                    if exist('eta', 'var')
                        total_targets = target_string(target, idx, k_target, total_targets, eta);     
                    else
                        total_targets = target_string(target, idx, k_target, total_targets);     
                    end 
                end 
                
                % every internal point has springs up one index in j and k 
                if is_internal(j,k,N)
                    
                    % nbr at j+1,k
                    rest_len = ref_frac * norm(R(:,j+1,k) - R(:,j,k)); 
                    kappa = k_rel / rest_len;         
                    nbr_idx = global_idx + vertex_index_offset(j+1,k,N);
                    total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    
                    % nbr at j,k+1
                    rest_len = ref_frac * norm(R(:,j,k+1) - R(:,j,k)); 
                    kappa = k_rel / rest_len;         
                    nbr_idx = global_idx + vertex_index_offset(j,k+1,N);
                    total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs);
                end 
                    
                pts_placed = pts_placed + 1; 
            end 
        end 
    end
    
    params.indices_global = indices_global; 
    global_idx = global_idx + pts_placed;

end 



function [global_idx, total_vertices, total_springs] = ...
                add_chordae_tree(params, spring, vertex, global_idx, total_vertices, total_springs, k_rel)

    % Adds chordae tree to IBAMR format files 
    % No targets here, so files and and count not included 
                        
    if ~isfield(params, 'chordae') || isempty(params.chordae)
        error('cannot place chordae with empty array'); 
    end 
    
    [X,alpha,beta,N,p_0,R,ref_frac,chordae] = unpack_params(params); 
    
    [m N_chordae] = size(chordae.C_left); 
    
    [C_left, C_right, left_papillary, right_papillary, Ref_l, Ref_r] = unpack_chordae(chordae); 
    
    
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
                
                % zero index means papillary muscles, which are zero and one by convention 
                if left_side 
                    nbr_idx = 0; 
                else 
                    nbr_idx = 1; 
                end
                
            else
                nbr_idx = parent + global_idx - 1; % subtract one, ibamr zero indexed  
                if ~left_side
                    nbr_idx = nbr_idx + N_chordae; 
                end 
            end 
            
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr R_nbr k_val j_nbr k_nbr] = get_nbr_chordae(params, i, parent, left_side); 
            
            rest_len = ref_frac * norm(R_nbr - Ref(:,i)); 
            kappa = k_rel * k_val / rest_len; 

            % list nbr index first because nbr is parent and has lower index
            total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs);        

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
            
            % just ignore the 
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
                            place_rays(params, filter_params, L, spring, vertex, target, ...
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

    r = filter_params.r; 
    h = filter_params.h; 
    N = params.N; 
    
    % max norm for included rays 
    % limit = L - (2*L / (2*N)); % leave points half mesh width from edge 
    
    if ~isfield(params, 'indices_global')
        error('Must place leaflets before placing rays'); 
    end 
    
    for k=1:N+1
        
        % working with the ring coordinates here 
        j = (N+2) - k; 
            
        pt_ring = params.X(:,j,k); 
        
        % only get a fiber if the previous point is included in the leaflet  
        neighbors = []; 
        if j > 1
            neighbors = [params.X(:,j-1,k), neighbors] ; 
        end 
        if k > 1
            neighbors = [params.X(:,j,k-1), neighbors] ; 
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
            nbr_idx = params.indices_global(j,k); 

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
    
    if abs(val(3)) > eps
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


function R = rotation_matrix_z(theta)
    % 
    % Rotation matrix counter clockwise by theta 
    % around z axis 

    R = [cos(theta) -sin(theta) 0; 
         sin(theta)  cos(theta) 0; 
         0           0          1]; 
end 


function R = rotation_matrix_y(theta)
    % 
    % Rotation matrix counter clockwise by theta 
    % around z axis 

    R = [cos(theta)  0 -sin(theta); 
         0           1  0         ; 
         sin(theta)  0  cos(theta)]; 
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



function [global_idx, total_vertices, total_lagrangian_placed] = place_lagrangian_tracers(global_idx, total_vertices, vertex, n_lagrangian_tracers, L, filter_params, double_z)
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
    
    z_extra_offset = filter_params.h; 
    
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





