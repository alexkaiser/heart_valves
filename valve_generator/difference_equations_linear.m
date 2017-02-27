function [F_leaflet F_chordae_left F_chordae_right] = difference_equations_linear(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Uses linear constitutive laws 
    % Requires rest lengths and spring constants to be previously
    % 
    % Input
    %     leaflet    Current parameters 
    %
    % Output
    %     F        Values of all difference equation, 3 by triangular array 
    % 

    X_current          = leaflet.X; 
    p_0                = leaflet.p_0; 
    C_left             = leaflet.chordae.C_left; 
    C_right            = leaflet.chordae.C_right; 
    chordae_idx_left   = leaflet.chordae_idx_left; 
    chordae_idx_right  = leaflet.chordae_idx_right;
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    du                 = leaflet.du; 
    dv                 = leaflet.dv; 
    is_internal        = leaflet.is_internal; 
    
    R_u = leaflet.R_u;
    k_u = leaflet.k_u;
    R_v = leaflet.R_v;
    k_v = leaflet.k_v;

    R_free_edge_left   = leaflet.R_free_edge_left;
    k_free_edge_left   = leaflet.k_free_edge_left;
    R_free_edge_right  = leaflet.R_free_edge_right;
    k_free_edge_right  = leaflet.k_free_edge_right;

    F_leaflet = zeros(size(X_current)); 

    [m N_chordae] = size(C_left); 


    free_edge_idx_left  = leaflet.free_edge_idx_left; 
    free_edge_idx_right = leaflet.free_edge_idx_right;

    for left_side = [true, false]

        if left_side
            free_edge_idx = free_edge_idx_left; 
            chordae_idx = chordae_idx_left; 
            C = C_left; 
            k_free_edge = k_free_edge_left; 
            R_free_edge = R_free_edge_left; 
        else 
            free_edge_idx = free_edge_idx_right; 
            chordae_idx = chordae_idx_right;
            C = C_right; 
            k_free_edge = k_free_edge_right; 
            R_free_edge = R_free_edge_right; 
        end 

        for i=1:size(free_edge_idx, 1)

            F_tmp = zeros(3,1);

            % left free edge has spring connections up and right on both leaflets
            j = free_edge_idx(i,1);
            k = free_edge_idx(i,2);

            X = X_current(:,j,k); 

            % interior neighbor is right in j on left side, 
            % left in j on right side 
            if left_side
                j_nbr = j + 1;
            else 
                j_nbr = j - 1;
            end 
            k_nbr = k;

            % Anterior circumferential 
            X_nbr = X_current(:,j_nbr,k_nbr); 
            
            tension = tension_linear(X,X_nbr,R_u(j,k),k_u(j,k)); 
            F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            % interior neighbor is up in k, always 
            j_nbr = j;     
            k_nbr = k+1; 

            % Anterior radial
            X_nbr = X_current(:,j_nbr,k_nbr); 
            tension = tension_linear(X,X_nbr,R_v(j,k),k_v(j,k)); 
            F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            % current node has a chordae connection
            if chordae_idx(j,k)

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k) + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = C(:,idx_chordae);
                
                tension = tension_linear(X,X_nbr,R_free_edge(i),k_free_edge(i));
                
                F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            else
                error('free edge point required to have chordae connection'); 
            end

            F_leaflet(:,j,k) = F_tmp; 

        end 

    end 


    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) && (~chordae_idx_left(j,k)) && (~chordae_idx_right(j,k))

                X = X_current(:,j,k); 

                F_tmp = zeros(3,1);

                % pressure term first  
                if p_0 ~= 0
                    F_tmp = F_tmp + (p_0 / (4*du*dv)) * cross(X_current(:,j+1,k) - X_current(:,j-1,k), X_current(:,j,k+1) - X_current(:,j,k-1));                     
                end 

                % u type fibers 
                for j_nbr = [j-1,j+1]

                    k_nbr = k; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 
                    
                    tension = tension_linear(X,X_nbr,R_u(j,k),k_u(j,k));                    
                    F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
                end 

                % v type fibers 
                for k_nbr = [k-1,k+1]
                    
                    j_nbr = j; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 
                    
                    tension = tension_linear(X,X_nbr,R_v(j,k),k_v(j,k)); 
                    F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
                end 

                F_leaflet(:,j,k) = F_tmp;

            end
        end 
    end
    

    % chordae internal terms 
    F_chordae_left  = zeros(size(C_left )); 
    F_chordae_right = zeros(size(C_right)); 

    [m N_chordae] = size(C_left); 

    for left_side = [true false];  

        if left_side
            C = C_left; 
        else 
            C = C_right; 
        end
        
        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference length and spring constants
                % routine handles unpacking and pulling correct constants 
                [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 
                
                tension = tension_linear(C(:,i),nbr,R_nbr,k_val); 
                
                tension_by_tangent = tension * (nbr - C(:,i)) / norm(nbr - C(:,i));  

                if left_side
                    F_chordae_left(:,i)  = F_chordae_left(:,i)  + tension_by_tangent; 
                else 
                    F_chordae_right(:,i) = F_chordae_right(:,i) + tension_by_tangent; 
                end 

            end 

        end 
    end 

end 



