function leaflet_linear = set_rest_lengths_and_constants_linear(leaflet, strain)
% 
% Assignes spring constants and rest lengths such that the current 
% valve configuration has uniform strain as specified here 
% 
% Tensions are computed from the configuration 
% 


X_current           = leaflet.X; 
alpha               = leaflet.alpha; 
beta                = leaflet.beta; 
C_left              = leaflet.chordae.C_left; 
C_right             = leaflet.chordae.C_right; 
k_0                 = leaflet.chordae.k_0; 
chordae_idx_left    = leaflet.chordae_idx_left; 
chordae_idx_right   = leaflet.chordae_idx_right;
j_max               = leaflet.j_max; 
k_max               = leaflet.k_max; 
du                  = leaflet.du; 
dv                  = leaflet.dv; 
is_internal         = leaflet.is_internal; 
free_edge_idx_left  = leaflet.free_edge_idx_left; 
free_edge_idx_right = leaflet.free_edge_idx_right;

[m N_chordae] = size(C_left); 

R_u = zeros(j_max, k_max); 
k_u = zeros(j_max, k_max); 
R_v = zeros(j_max, k_max); 
k_v = zeros(j_max, k_max); 

n_free_edge_left   = size(free_edge_idx_left, 1); 
n_free_edge_right  = size(free_edge_idx_left, 1); 

R_free_edge_left   = zeros(n_free_edge_left,  1); 
k_free_edge_left   = zeros(n_free_edge_left,  1); 
R_free_edge_right  = zeros(n_free_edge_right, 1); 
k_free_edge_right  = zeros(n_free_edge_right, 1);

R_chordae_left  = zeros(N_chordae,1); 
k_chordae_left  = zeros(N_chordae,1); 
R_chordae_right = zeros(N_chordae,1); 
k_chordae_right = zeros(N_chordae,1); 


% repulsive potential coefficients, if used 
if isfield(leaflet, 'repulsive_potential') && leaflet.repulsive_potential
    repulsive_potential         = true; 
    power                       = leaflet.repulsive_power; 
    c_repulsive_circumferential = leaflet.c_repulsive_circumferential; 
    c_repulsive_radial          = leaflet.c_repulsive_radial; 
    c_repulsive_chordae         = leaflet.c_repulsive_chordae; 
else 
    repulsive_potential         = false; 
    power                       = 1; 
    c_repulsive_circumferential = 0.0; 
    c_repulsive_radial          = 0.0; 
    c_repulsive_chordae         = 0.0; 
end 

if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
    decreasing_tension = true; 
    c_dec_tension_circumferential = leaflet.c_dec_tension_circumferential; 
    c_dec_tension_radial          = leaflet.c_dec_tension_radial; 
    c_dec_tension_chordae         = leaflet.c_dec_tension_chordae; 
else 
    decreasing_tension = false; 
    c_dec_tension_circumferential = 0.0; 
    c_dec_tension_radial          = 0.0; 
    c_dec_tension_chordae         = 0.0; 
end 
    

[m N_chordae] = size(C_left); 



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
            
        % spring and rest length indices are always at the minimum value  
        j_spr = min(j, j_nbr); 
        k_spr = min(k, k_nbr); 
            
        % Anterior circumferential 
        X_nbr = X_current(:,j_nbr,k_nbr); 

        % Multiply tension by dv to get a force,
        % rather than a force density, here 
        tension = alpha * dv; 

        if repulsive_potential
            tension = tension - alpha * dv * c_repulsive_circumferential * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
        end 

        if decreasing_tension
            tension = tension + alpha * dv * tension_decreasing(X, X_nbr, du, c_dec_tension_circumferential) ; 
        end 

        [k_u(j_spr,k_spr) R_u(j_spr,k_spr)] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 

        % interior neighbor is up in k, always 
        j_nbr = j;     
        k_nbr = k+1; 
        j_spr = min(j, j_nbr); 
        k_spr = min(k, k_nbr); 

        % Anterior radial
        X_nbr = X_current(:,j_nbr,k_nbr); 
        tension = beta * du; 

        if repulsive_potential
            tension = tension - beta * du * c_repulsive_radial * dv^2 * power * 1/norm(X_nbr-X)^(power+1); 
        end 

        if decreasing_tension
            tension = tension + beta * du * tension_decreasing(X, X_nbr, dv, c_dec_tension_radial) ; 
        end

        [k_v(j_spr,k_spr) R_v(j_spr,k_spr)] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 

        % current node has a chordae connection
        if chordae_idx(j,k)

            kappa = k_0;

            % index that free edge would have if on tree
            % remember that leaves are only in the leaflet
            leaf_idx = chordae_idx(j,k) + N_chordae;

            % then take the parent index of that number in chordae variables
            idx_chordae = floor(leaf_idx/2);

            X_nbr = C(:,idx_chordae);
            tension = kappa; 

            if repulsive_potential
                tension = tension - kappa * c_repulsive_chordae * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
            end 

            if decreasing_tension
                tension = tension + kappa * tension_decreasing(X, X_nbr, du, c_dec_tension_chordae) ; 
            end

            if left_side 
                [k_free_edge_left(i), R_free_edge_left(i)]   = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 
            else 
                [k_free_edge_right(i), R_free_edge_right(i)] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 
            end 
            
        else
            error('free edge point required to have chordae connection'); 
        end
    end 
    
end 


% Internal leaflet part 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) && (~chordae_idx_left(j,k)) && (~chordae_idx_right(j,k))

            X = X_current(:,j,k); 

            % u type fibers 
            % set constants in up direction only here 
            for j_nbr = [j-1,j+1] 
                
                k_nbr = k; 

                j_spr = min(j, j_nbr); 
                k_spr = min(k, k_nbr);
                
                X_nbr = X_current(:,j_nbr,k_nbr); 

                tension = alpha; 

                if repulsive_potential
                    tension = tension - alpha * c_repulsive_circumferential * du^2 * power * 1/norm(X_nbr-X)^(power+1); 
                end 

                if decreasing_tension
                    tension = tension + alpha * tension_decreasing(X, X_nbr, du, c_dec_tension_circumferential) ; 
                end 

                % Here tension is a force per unit length 
                % So we must multiply by a legnth element to get force 
                % Take the opposing length element 
                tension = dv * tension; 

                [k_u(j_spr,k_spr) R_u(j_spr,k_spr)] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 
                
            end 


            % v type fibers 
            for k_nbr = [k-1,k+1] 
                
                j_nbr = j; 
                
                j_spr = min(j, j_nbr); 
                k_spr = min(k, k_nbr);
                
                X_nbr = X_current(:,j_nbr,k_nbr); 

                tension = beta; 

                if repulsive_potential
                    tension = tension - beta * c_repulsive_radial * dv^2 * power * 1/norm(X_nbr-X)^(power+1); 
                end 

                if decreasing_tension
                    tension = tension + beta * tension_decreasing(X, X_nbr, dv, c_dec_tension_radial) ; 
                end

                tension = du * tension; 

                [k_v(j_spr,k_spr) R_v(j_spr,k_spr)] = get_rest_len_and_spring_constant_linear(X, X_nbr, tension, strain); 
                
            end 
        end
    end 
end
 


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

%         for nbr_idx = [left,right,parent] 
         for nbr_idx = [parent] 
            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, left_side); 

            tension = k_val; 

            if repulsive_potential
                tension = tension - k_val * c_repulsive_chordae * du^2 * power * 1/norm(nbr - C(:,i))^(power+1); 
            end 

            if decreasing_tension
                tension = tension + k_val * tension_decreasing(C(:,i), nbr, du, c_dec_tension_chordae) ; 
            end

            if left_side
                [k_chordae_left(i)  R_chordae_left(i) ] = get_rest_len_and_spring_constant_linear(C(:,i), nbr, tension, strain); 
            else 
                [k_chordae_right(i) R_chordae_right(i)] = get_rest_len_and_spring_constant_linear(C(:,i), nbr, tension, strain); 
            end 

        end          
    end 
end 

% Copy all basic data structures 
leaflet_linear = leaflet; 

% Add new information 
leaflet_linear.R_u = R_u;
leaflet_linear.k_u = k_u;
leaflet_linear.R_v = R_v;
leaflet_linear.k_v = k_v;

leaflet_linear.R_free_edge_left   = R_free_edge_left;
leaflet_linear.k_free_edge_left   = k_free_edge_left;
leaflet_linear.R_free_edge_right  = R_free_edge_right;
leaflet_linear.k_free_edge_right  = k_free_edge_right;

leaflet_linear.chordae.Ref_l  = R_chordae_left;
leaflet_linear.chordae.k_l    = k_chordae_left;
leaflet_linear.chordae.Ref_r  = R_chordae_right;
leaflet_linear.chordae.k_r    = k_chordae_right;

leaflet_linear.diff_eqns = @difference_equations_linear; 
leaflet_linear.jacobian = @build_jacobian_linear; 



