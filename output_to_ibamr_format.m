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
    k_target = k_rel; 

    % output the left and right papillary as the first two vertices and targets 
    left_papillary  = [0; -filter_params_posterior.a; 0]; 
    total_vertices = vertex_string(vertex, left_papillary, total_vertices); 
    total_targets = target_string(target, global_idx, k_target, total_targets);     
    global_idx = global_idx + 1; 
    
    right_papillary = [0;  filter_params_posterior.a; 0]; 
    total_vertices = vertex_string(vertex, right_papillary, total_vertices); 
    total_targets = target_string(target, global_idx, k_target, total_targets);     
    global_idx = global_idx + 1;     
        
    
    [X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params_posterior); 
    
    pts_placed = 0; 
    
    for k=1:N+1
        for j=1:N+1
            if (j+k) <= (N+2)
                
                idx = global_idx + vertex_index_offset(j,k,N); 
             
                j
                k
                idx 
                
                
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
                    'j+1 nbr'
                    nbr_idx = global_idx + vertex_index_offset(j+1,k,N)
                    total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs); 
                    
                    % nbr at j,k+1
                    rest_len = ref_frac * norm(R(:,j,k+1) - R(:,j,k)); 
                    kappa = k_rel / rest_len;         
                    nbr_idx = global_idx + vertex_index_offset(j,k+1,N)
                    'k+1 nbr'
                    total_springs = spring_string(spring, idx, nbr_idx, kappa, rest_len, total_springs);
                end 
                    
                pts_placed = pts_placed + 1; 
            end 
        end 
    end 
    
    global_idx = global_idx + pts_placed; 
    



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
    movefile(strcat(file_name, '.tmp'),file_name); 

end 





















