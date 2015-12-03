function [] = output_to_ibamr_format(base_name, L, ratio, params_posterior, filter_params_posterior, params_anterior, filter_params_anterior)
    % 
    % Outputs the current configuration of the leaflets to IBAMR format
    % Spring constants are computed in dimensional form 
    % 
    %
    % Input: 
    %    base_name                  File base name
    %    L                          Outputs extra mesh to use a [-L,L]^3 cube (not implemented)
    %    ratio                      Ratio of pressure to nondimensionalized spring constant   
    %    params_posterior           Parameters for various leaflets 
    %    filter_params_posterior
    %    params_anterior
    %    filter_params_anterior 
    % 
    % Output: 
    %    Files written in IBAMR format 
    % 

    vertex = fopen(strcat(base_name, '.vertex'), 'w'); 
    spring = fopen(strcat(base_name, '.spring'), 'w'); 
    target = fopen(strcat(base_name, '.target'), 'w'); 

    % keep one global index through the whole thing 
    global_idx = 0;
    
    % just count the number of vertices and strings throughout 
    total_vertices = 0; 
    total_springs  = 0; 
    total_targets  = 0; 

    % compute some needed constants 
    MMHG_TO_CGS = 1333.22368; 
    p_cgs = 10 * MMHG_TO_CGS; 

    % value of the relative spring constant is determined by the ratio 
    k_rel = p_cgs / ratio; 
    
    % turn these up because they are supposed to really be boundary conditions 
    % but we are using a penalty method here 
    k_target = 100 * k_rel; 

    % output the left and right papillary as the first two vertices and targets 
    left_papillary  = [0; -filter_params_posterior.a; 0]; 
    total_vertices = vertex_string(vertex, left_papillary, total_vertices); 
    total_targets = target_string(target, global_idx, k_target, total_targets);     
    global_idx = global_idx + 1; 
    
    right_papillary = [0;  filter_params_posterior.a; 0]; 
    total_vertices = vertex_string(vertex, right_papillary, total_vertices); 
    total_targets = target_string(target, global_idx, k_target, total_targets);     
    global_idx = global_idx + 1;     
        
    % leaflets 
    [global_idx, total_vertices, total_springs, total_targets] = ...
        add_leaflet(params_posterior, left_papillary, right_papillary, spring, vertex, target, ...
                    global_idx, total_vertices, total_springs, total_targets, k_rel, k_target); 
        
    [global_idx, total_vertices, total_springs, total_targets] = ...
         add_leaflet(params_anterior, left_papillary, right_papillary, spring, vertex, target, ...
                     global_idx, total_vertices, total_springs, total_targets, k_rel, k_target);    

    % flat part of mesh 
    L = 3.0; 
    r = filter_params_posterior.r; 
    N_ring = 2 * params_anterior.N;
    h = filter_params_posterior.h; 
    ref_frac = params_posterior.ref_frac; 
    
    [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_net(r, h, L, N_ring, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, ref_frac); 
                        

    % clean up files with totals 
    fclose(vertex); 
    fclose(spring); 
    fclose(target); 

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

function total_targets = target_string(target, idx, kappa, total_targets)
    % prints a target format string to target file 
    fprintf(target, '%d\t %.14f\n', idx, kappa); 
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

function [global_idx, total_vertices, total_springs, total_targets] = ...
                add_leaflet(params, left_papillary, right_papillary, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target)
                                                                                
% Places all main data into IBAMR format for this leaflet
% Updates running totals on the way 
                                                                                
                                                                                

    [X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params); 
    
    pts_placed = 0; 
    
    for k=1:N+1
        for j=1:N+1
            if (j+k) <= (N+2)
                
                idx = global_idx + vertex_index_offset(j,k,N);    
                
                total_vertices = vertex_string(vertex, X(:,j,k), total_vertices); 
                
                % if j==1, connect to the left papillary 
                if j==1
                    nbr_idx = 0; 
                    rest_len = ref_frac * norm(left_papillary - R(:,j,k)); 
                    kappa = k_rel / rest_len; 
                    
                    total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs); 
                end 
                
                % if k==1, connect to the right papillary 
                if k==1
                    nbr_idx = 1; 
                    rest_len = ref_frac * norm(right_papillary - R(:,j,k)); 
                    kappa = k_rel / rest_len; 
                    
                    total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs); 
                end                 
                
                % if on boundary, this is a target point 
                if (j+k) == (N+2)
                    total_targets = target_string(target, idx, k_target, total_targets);     
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
    
    global_idx = global_idx + pts_placed; 

end 


function [global_idx, total_vertices, total_springs, total_targets] = ...
                            place_net(r, h, L, N, spring, vertex, target, ...
                            global_idx, total_vertices, total_springs, total_targets, k_rel, k_target, ref_frac)
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
    %     global_idx       Running totals  
    %     total_vertices 
    %     total_springs 
    %     total_targets 
    % 

    % mesh spacing on the valve ring 
    ds = 2*pi / N; 

    % maximum number of points in radial direction 
    % this gets a sqrt(2) because we are placing the net in a square 
    M = floor(sqrt(2) * L / ds); 

    points = zeros(3,N,M);

    % just keep a list of valid indices, mark NAN if out of physical bounds 
    indices = zeros(N,M); 
    idx = 0; 
    points_placed = 0; 
    

    rad = r; 
    
    for k=1:M
        for j=1:N
            
            theta = (j-1) * ds; 
            coords = [rad*cos(theta); rad*sin(theta); h]; 
            
            % if one norm is less than L, then the point is within the domain  
            if norm(coords, inf) < L 
                points(:,j,k) = coords; 
                indices(j,k) = idx; 
                idx = idx + 1;                 
            else 
                points(:,j,k) = NaN * ones(3,1);
                indices(j,k) = NaN; 
            end     
            
        end 
        
        rad = rad + ds; 
    end 

    
    % below the first possible point 
    idx = -1; 
    
    for k=1:M
        for j=1:N
            
            % just ignore the 
            if ~isnan(indices(j,k))
 
                last_idx = idx; 
                idx = indices(j,k) + global_idx; 
                
                if last_idx >= idx
                    error('should always be placing points in order, something wrong'); 
                end 
                
                % every valid vertex is a target here 
                total_vertices = vertex_string(vertex, points(:,j,k), total_vertices); 
                points_placed = points_placed + 1; 
                total_targets = target_string(target, idx, k_target, total_targets);     
                
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
                if (k+1) < M
                    if ~isnan(indices(j,k+1))
                        rest_len = ref_frac * norm(points(:,j,k) - points(:,j,k+1)); 
                        kappa = k_rel / rest_len; 
                        nbr_idx = indices(j,k+1) + global_idx; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    end 
                end 
                
                % no periodic direction in radial direction (k)
                
            end 
           
        end 
    end 
    
    global_idx = global_idx + points_placed; 

end 









