function [] = output_leaflet_to_xyz_format(leaflet, springs_needed)
    % 
    % Outputs the current configuration of the leaflets to IBAMR format
    % Spring constants are computed in dimensional form 
    % 
    %
    % Input: 
    %    base_name                  File base name
    %    L                          Outputs extra mesh to use a [-L,L]^3 cube
    %    ratio                      Ratio of pressure to nondimensionalized spring constant   
    %    params_leaflet           Parameters for various leaflets 
    %    filter_params_leaflet
    %    params_leaflet
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
    
    base_name       = leaflet.base_name; 
    
    vertex_name = sprintf('%s_%.10d.xyz', base_name, leaflet.frame); 
    spring_name = strcat(base_name, '.spring'); 
    
    vertex = fopen(vertex_name, 'w'); 
    fprintf(vertex, '\n'); % always a wasted line at beginning 
    
    if springs_needed
        spring = fopen(spring_name, 'w'); 
    else
        spring = false;  
    end
    
    % keep one global index through the whole thing 
    global_idx = 0;
    
    % just count the number of vertices and strings throughout 
    total_vertices = 0; 
    total_springs  = 0; 
  
    k_rel_leaflet = 1.0; 
    

    if leaflet.reflect_x
        leaflet.X(1,:,:)                   = -leaflet.X(1,:,:); 
        leaflet.R(1,:,:)                   = -leaflet.R(1,:,:); 
        leaflet.left_papillary(1)          = -leaflet.left_papillary(1); 
        leaflet.right_papillary(1)         = -leaflet.right_papillary(1); 
        leaflet.chordae.C_left (1,:,:)     = -leaflet.chordae.C_left (1,:,:);
        leaflet.chordae.C_right(1,:,:)     = -leaflet.chordae.C_right(1,:,:);
        leaflet.chordae.left_papillary(1)  = -leaflet.chordae.left_papillary(1); 
        leaflet.chordae.right_papillary(1) = -leaflet.chordae.right_papillary(1); 
    end 
    
    
    % check for consistency in chordae , all data structures must match 
    if ~(    all(leaflet.left_papillary   == leaflet.chordae.left_papillary)   ...
          && all(leaflet.right_papillary  == leaflet.chordae.right_papillary))
        error('chordae are inconsistent'); 
    end 
    
    % always output current setup here 
    leaflet.R               = leaflet.X; 
    leaflet.chordae.Ref_l   = leaflet.chordae.C_left;  
    leaflet.chordae.Ref_r   = leaflet.chordae.C_right;  
    leaflet.ref_frac        = 1.0; 
    
    
    % output the left and right papillary as the first two vertices and targets
    left_papillary  = leaflet.left_papillary; 
    leaflet.left_papillary_idx = global_idx; 
    total_vertices  = vertex_string(vertex, left_papillary, total_vertices); 
    global_idx      = global_idx + 1; 

    right_papillary = leaflet.right_papillary; 
    leaflet.right_papillary_idx = global_idx; 
    total_vertices  = vertex_string(vertex, right_papillary, total_vertices); 
    global_idx      = global_idx + 1;     

    % leaflet 
    [global_idx, total_vertices, total_springs, leaflet] = ...
        add_leaflet(leaflet, spring, vertex, ...
                     global_idx, total_vertices, total_springs, k_rel_leaflet);   

    % leaflet chordae 
    leaflet.chordae.left_papillary_idx  = leaflet.left_papillary_idx; 
    leaflet.chordae.right_papillary_idx = leaflet.right_papillary_idx; 

    [global_idx, total_vertices, total_springs] = ...
            add_chordae_tree(leaflet, spring, vertex, global_idx, total_vertices, total_springs, k_rel_leaflet);  
                

    % clean up files with totals 
    fclose(vertex); 
    
    prepend_line_with_int(vertex_name, total_vertices); 
    
    if springs_needed 
        fclose(spring); 
        prepend_line_with_int(spring_name, total_springs);
    end 
end 



% nest this function so it can access the z increment
% this is bad practice and should be removed 
function total_vertices = vertex_string(vertex, coords, total_vertices)
    % prints formatted string for current vertex to vertex file   
    
    fprintf(vertex, '? %.14f\t %.14f\t %.14f\n', coords(1), coords(2), coords(3)); 
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


function [global_idx, total_vertices, total_springs, leaflet] = ...
                add_leaflet(leaflet, spring, vertex, ...
                            global_idx, total_vertices, total_springs, k_rel)
                      
    % Places all main data into IBAMR format for this leaflet
    % Updates running totals on the way 

    % Unpack needed data 
    X                 = leaflet.X; 
    R                 = leaflet.R; 
    j_max             = leaflet.j_max; 
    k_max             = leaflet.k_max; 
    is_internal       = leaflet.is_internal;
    point_idx_with_bc = leaflet.point_idx_with_bc; 
    is_bc             = leaflet.is_bc; 
    chordae           = leaflet.chordae;
    chordae_idx_left  = leaflet.chordae_idx_left; 
    chordae_idx_right = leaflet.chordae_idx_right; 
    ref_frac          = leaflet.ref_frac; 
    alpha             = leaflet.alpha; 
    beta              = leaflet.beta; 
    
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
                
                if spring 
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

                        kappa = k_rel * chordae.k_0 / rest_len; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                        
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

                        kappa = k_rel * chordae.k_0 / rest_len; 
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 

                    end 


                    % springs in leaflet, only go in up direction 
                    j_nbr = j + 1; 
                    k_nbr = k; 
                    if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))

                        rest_len = ref_frac * norm(R(:,j_nbr,k_nbr) - R(:,j,k)); 
                        nbr_idx = global_idx + point_idx_with_bc(j_nbr,k_nbr);

                        kappa = alpha * k_rel / rest_len;         
                        total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 

                    end 

                    % springs in leaflet, only go in up direction 
                    j_nbr = j; 
                    k_nbr = k + 1; 
                    if (j_nbr <= j_max) && (k_nbr <= k_max) && (is_internal(j_nbr,k_nbr) || is_bc(j_nbr, k_nbr))

                        rest_len = ref_frac * norm(R(:,j_nbr,k_nbr) - R(:,j,k)); 
                        nbr_idx = global_idx + point_idx_with_bc(j_nbr,k_nbr);

                        kappa = beta * k_rel / rest_len;         
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
                add_chordae_tree(leaflet, spring, vertex, global_idx, total_vertices, total_springs, k_rel)

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
            
            if spring 
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

                kappa = k_rel * k_val / rest_len; 
                % list nbr index first because nbr is parent and has lower index
                total_springs = spring_string(spring, nbr_idx, idx, kappa, rest_len, total_springs);        

            end 
            
        end 
        
    end 
    
    global_idx = global_idx + 2*N_chordae; 

end 
                        
                        










