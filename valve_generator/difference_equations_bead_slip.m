function [F_leaflet F_chordae_left F_chordae_right] = difference_equations_bead_slip(leaflet)
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
    R_current          = leaflet.R; 
    p_0                = leaflet.p_0; 
    alpha              = leaflet.alpha; 
    beta               = leaflet.beta; 
    ref_frac           = leaflet.ref_frac; 
    C_left             = leaflet.chordae.C_left; 
    C_right            = leaflet.chordae.C_right; 
    Ref_l              = leaflet.chordae.Ref_l; 
    Ref_r              = leaflet.chordae.Ref_r; 
    k_0                = leaflet.chordae.k_0; 
    chordae_idx_left   = leaflet.chordae_idx_left; 
    chordae_idx_right  = leaflet.chordae_idx_right;
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    du                 = leaflet.du; 
    dv                 = leaflet.dv; 
    is_internal        = leaflet.is_internal; 


    F_leaflet = zeros(size(X_current)); 


    [m N_chordae] = size(C_left); 


    S_left  = zeros(k_max-1,1); 
    S_right = zeros(k_max-1,1);     
    T       = zeros(j_max,1); 


    free_edge_idx_left  = leaflet.free_edge_idx_left; 
    free_edge_idx_right = leaflet.free_edge_idx_right;

    for left_side = [true, false]

        if left_side
            free_edge_idx = free_edge_idx_left; 
            chordae_idx = chordae_idx_left; 
            C = C_left; 
            Ref = Ref_l;
        else 
            free_edge_idx = free_edge_idx_right; 
            chordae_idx = chordae_idx_right;
            C = C_right; 
            Ref = Ref_r;
        end 

        for i=1:size(free_edge_idx, 1)

            F_tmp = zeros(3,1);

            % left free edge has spring connections up and right on both leaflets
            j = free_edge_idx(i,1);
            k = free_edge_idx(i,2);

            X = X_current(:,j,k); 
            R = R_current(:,j,k);

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
            R_nbr = R_current(:,j_nbr,k_nbr); 

            if left_side
                S_left(k)  = tension_linear(X, X_nbr, R, R_nbr, alpha, ref_frac); 
                F_tmp = F_tmp + S_left(k) * (X_nbr-X)/norm(X_nbr-X); 
            else
                S_right(k) = tension_linear(X, X_nbr, R, R_nbr, alpha, ref_frac); 
                F_tmp = F_tmp + S_right(k) * (X_nbr-X)/norm(X_nbr-X);             
            end

            % interior neighbor is up in k, always 
            j_nbr = j;     
            k_nbr = k+1; 

            % Anterior radial
            X_nbr = X_current(:,j_nbr,k_nbr); 
            R_nbr = R_current(:,j_nbr,k_nbr); 
            T(j) = tension_linear(X, X_nbr, R, R_nbr, beta, ref_frac); 
            F_tmp = F_tmp + T(j) * (X_nbr-X)/norm(X_nbr-X); 

            % current node has a chordae connection
            if chordae_idx(j,k)

                kappa = k_0;

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k) + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = C(:,idx_chordae);
                R_nbr = Ref(:,idx_chordae);

                F_tmp = F_tmp + tension_linear(X,X_nbr,R,R_nbr,kappa,ref_frac) * (X_nbr-X)/norm(X_nbr-X); 

            else
                error('free edge point required to have chordae connection'); 
            end

            F_leaflet(:,j,k) = F_tmp; 

        end 

    end 

    % Tensions are equalized from left to right 
    S = (S_left + S_right)/2.0;
 
    % Convert from units of force to force densities 
    S = S/dv;
    T = T/du;

    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) && (~chordae_idx_left(j,k)) && (~chordae_idx_right(j,k))

                X = X_current(:,j,k); 

                F_tmp = zeros(3,1);

                % pressure term first  
                if p_0 ~= 0
                    F_tmp = F_tmp + (p_0 / (du*dv)) * cross(X(:,j+1,k) - X(:,j-1,k), X(:,j,k+1) - X(:,j,k-1));                     
                end 

                if false
                % u type fibers 
                for j_nbr = [j-1,j+1]

                    k_nbr = k; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    F_tmp = F_tmp + S(k)/du * (X_nbr-X)/norm(X_nbr-X); 
                    % F_tmp = F_tmp + 1/du * (X_nbr-X)/norm(X_nbr-X); 

                end 
                end 

                if true
                % v type fibers 
                for k_nbr = [k-1,k+1]

                    j_nbr = j; 
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    F_tmp = F_tmp + T(j)/dv * (X_nbr-X)/norm(X_nbr-X); 

                end 
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
            Ref = Ref_l; 
        else 
            C = C_right; 
            Ref = Ref_r; 
        end

        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference coordinate and spring constants
                [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 

                tension = tension_linear_over_norm(C(:,i), nbr, Ref(:,i), R_nbr, k_val, ref_frac) * (nbr - C(:,i));  

                if left_side
                    F_chordae_left(:,i)  = F_chordae_left(:,i)  + tension; 
                else 
                    F_chordae_right(:,i) = F_chordae_right(:,i) + tension; 
                end 

            end 

        end 
    end 

end 



