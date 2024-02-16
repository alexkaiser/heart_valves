function [free_edge_length_single_loaded, ...
          free_edge_length_single_rest, ...
          portion_of_current_edge, ...
          portion_of_restlen_edge, ...
          strains] ... 
           = get_circ_edge_lengths(leaflet, N_each, k, X, R_u, debug_lengths)

    free_edge_length_single_loaded = 0; 
    free_edge_length_single_rest = 0; 
    
    portion_of_current_edge = zeros(N_each,1);
    portion_of_restlen_edge = zeros(N_each,1);
    
    if ~exist('debug_lengths', 'var')
        debug_lengths = false; 
    end 
    
    loaded_lens = zeros(N_each,1); 
    rest_lens = zeros(N_each,1); 
    
    for j=1:N_each
        % k passed in 

        j_nbr_tmp = j-1; 
        k_nbr_tmp = k; 
        [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = X(:,j,k);
        X_nbr = X(:,j_nbr,k_nbr); 

        loaded_lens(j) = norm(X_temp - X_nbr); 
        rest_lens(j) = R_u(j_spr,k_spr); 

        free_edge_length_single_loaded = free_edge_length_single_loaded + norm(X_temp - X_nbr);        
        free_edge_length_single_rest = free_edge_length_single_rest + R_u(j_spr,k_spr); 
        
        portion_of_current_edge(j) = free_edge_length_single_loaded;
        portion_of_restlen_edge(j) = free_edge_length_single_rest; 
    end
    
    strains = (loaded_lens ./ rest_lens) - 1; 
    if debug_lengths 
%         fprintf('in get_circ_edge_lengths')
%         loaded_lens 
%         rest_lens 
        strains 
%         fprintf('end debug output get_circ_edge_lengths')
    end 
    
    portion_of_current_edge = portion_of_current_edge / free_edge_length_single_loaded;
    portion_of_restlen_edge = portion_of_restlen_edge / free_edge_length_single_rest;
    
end 