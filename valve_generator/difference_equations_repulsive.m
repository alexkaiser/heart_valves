function [F_leaflet F_chordae_left F_chordae_right] = difference_equations_repulsive(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Requires reference configuration R 
    % 
    % Input
    %     leaflet    Current parameters 
    %
    % Output
    %     F        Values of all difference equation, 3 by triangular array 
    % 

    X_current          = leaflet.X; 
    alpha              = leaflet.alpha; 
    beta               = leaflet.beta;
    C_left             = leaflet.chordae.C_left; 
    C_right            = leaflet.chordae.C_right;
    k_0                = leaflet.chordae.k_0; 
    chordae_idx_left   = leaflet.chordae_idx_left; 
    chordae_idx_right  = leaflet.chordae_idx_right;
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    du                 = leaflet.du; 
    dv                 = leaflet.dv;
    is_internal        = leaflet.is_internal; 

    if ~leaflet.repulsive_potential
        error('Calling repulsive difference equations with flags off'); 
    end 
    
    p                  = leaflet.repulsive_power; 
    coeff              = leaflet.repulsive_coeff; 
    

    F_leaflet = zeros(size(X_current)); 

    [m N_chordae] = size(C_left); 


    free_edge_idx_left  = leaflet.free_edge_idx_left; 
    free_edge_idx_right = leaflet.free_edge_idx_right;

    for left_side = [true, false]

        if left_side
            free_edge_idx = free_edge_idx_left; 
            chordae_idx = chordae_idx_left; 
            C = C_left; 
        else 
            free_edge_idx = free_edge_idx_right; 
            chordae_idx = chordae_idx_right;
            C = C_right; 
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

            F_tmp = F_tmp - alpha * dv * coeff * p * (X_nbr-X)/norm(X_nbr-X)^(p+2); 

            % interior neighbor is up in k, always 
            j_nbr = j;     
            k_nbr = k+1; 

            % Anterior radial
            X_nbr = X_current(:,j_nbr,k_nbr); 
            F_tmp = F_tmp - beta * du * coeff * p * (X_nbr-X)/norm(X_nbr-X)^(p+2); 

            % current node has a chordae connection
            if chordae_idx(j,k)
                
                kappa = k_0;

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k) + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = C(:,idx_chordae);

                F_tmp = F_tmp - kappa * coeff * p * (X_nbr-X)/norm(X_nbr-X)^(p+2); 

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

                % u type fibers 
                for j_nbr = [j-1,j+1]

                    k_nbr = k; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    F_tmp = F_tmp - alpha/du * coeff * (du^2) * p * (X_nbr-X)/norm(X_nbr-X)^(p+2); 

                end 

                % v type fibers 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    F_tmp = F_tmp - beta/dv * coeff * (dv^2) * p * (X_nbr-X)/norm(X_nbr-X)^(p+2); 

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

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 

                tension = - k_val * coeff * p * (nbr - C(:,i)) / norm(nbr - C(:,i))^(p+2);  

                if left_side
                    F_chordae_left(:,i)  = F_chordae_left(:,i)  + tension; 
                else 
                    F_chordae_right(:,i) = F_chordae_right(:,i) + tension; 
                end 

            end 

        end 
    end 

end 



