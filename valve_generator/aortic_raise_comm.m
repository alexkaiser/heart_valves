function valve_with_reference = aortic_raise_comm(valve_with_reference)

    % dialation with preserved arc length 

    if isfield(valve_with_reference, 'comm_raise_normal_height')
        comm_raise_normal_height = valve_with_reference.comm_raise_normal_height; 
    else 
        error('must provide comm_raise_h1 if valve.raise_comm is true'); 
    end 

    leaflet = valve_with_reference.leaflets(1); 

    hc = valve_with_reference.skeleton.normal_height - valve_with_reference.skeleton.height_min_comm; 
    
    height_min_comm_new = comm_raise_normal_height - hc; 
    
    % compute annular length initial 
    [~, len_annulus_min_initial, len_annulus_each] = build_initial_fibers_aortic(leaflet, valve_with_reference); 
    
    power_initial_guess = 3; 
    
    options = optimset('Display','off','TolFun',1e-16);
    
    % this has two return values 
    len_annulus_tmp = @(power) build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm_new, power); 

    % hack to return only the next one 
    len_annulus = @(power) Out2(len_annulus_tmp, power); 
    
    len_annulus_minus_initial = @(power) abs(len_annulus(power) - len_annulus_min_initial); 
    
    use_fsolve = false; 
    use_fmincon = true; 
    
    if use_fsolve 
        [power_new, fval, exitflag, output] = fsolve(len_annulus_minus_initial,power_initial_guess,options);     
    elseif use_fmincon
        min_bd = 1; 
        [power_new, fval, exitflag, output] = fmincon(len_annulus_minus_initial,power_initial_guess,[],[],[],[],min_bd,[],[],options);     
    end 
        
    if power_new < 1
        warning('failed to fine power greater than 1 for interpolation')
    end 
    
    if exitflag < 1
        warning('Errors in optimization: power_new = %f, exitflag = %d, output = %s', power_new, exitflag, output.message); 
    end 
    
    % reset skeleton and fibers post optimization 
    % total height 
    valve_with_reference.skeleton.normal_height = height_min_comm_new + hc; 

    % height to bottom of comm 
    valve_with_reference.skeleton.height_min_comm = height_min_comm_new; 
        
    valve_with_reference.skeleton.power = power_new; 
    
    [X, len_annulus_new] = build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm_new, power_new); 
    
    tol = 1e-8;
    if norm(len_annulus_min_initial - len_annulus_new) > tol 
        warning('new annulus len not equal to previous after raising comm')
    end 
    
    redistribute_annulus_points = false; 
    if redistribute_annulus_points 
        % respace free edge points 
        N_each = leaflet.N_each; 
        j_max = leaflet.j_max; 
        debug_lengths = false; 
        k = 1; 
        [~, ~, portion_of_current_edge, portion_of_rest_edge] = get_circ_edge_lengths(leaflet, N_each, k, X, len_annulus_each, debug_lengths);

        X_annular_leaflet_1_with_wrap = [X(:,j_max,k), X(:, 1:N_each,k)]; 

        % spacing of points as fraction of arc length 
        interp_idx_annular = [0; portion_of_current_edge];  

        % interpolate as fraction of rest length 
        annular_edge_interp_points_respaced = interp1(interp_idx_annular, X_annular_leaflet_1_with_wrap', portion_of_rest_edge(1:N_each-1))'; 

        X(:,(1:N_each-1),k) = annular_edge_interp_points_respaced; 

        % other leaflet by rotation 
        X(:, (N_each+1:j_max-1),k) = rotation_matrix_z(pi) * annular_edge_interp_points_respaced;     
        
        debug_redistribute = false; 
        if debug_redistribute
            % minimum ring 
            len_annulus_each_updated = zeros(j_max,1); 
            k = 1; 
            for j = 1:j_max

                [valid, j_nbr, k_nbr] = get_indices(leaflet, j, k, j+1, k); 
                if ~valid
                    error('failed to find valid index'); 
                end 
                len_annulus_each_updated(j) = norm(X(:,j_nbr,k_nbr) - X(:,j,k)); 
            end 
                        
            len_annulus_each
            len_annulus_each_updated
        end 
        
        % minimum ring 
        len_annulus_new_reinterp = 0; 
        k = 1; 
        for j = 1:j_max

            [valid, j_nbr, k_nbr] = get_indices(leaflet, j, k, j+1, k); 
            if ~valid
                error('failed to find valid index'); 
            end 

            len_annulus_new_reinterp = len_annulus_new_reinterp + norm(X(:,j_nbr,k_nbr) - X(:,j,k));     
        end 
        
        if norm(len_annulus_min_initial - len_annulus_new_reinterp) > tol 
            warning('new annulus len not equal to previous after raising comm, initial = %f, reinterp = %f, diff = %e\n', ...
                len_annulus_min_initial, len_annulus_new_reinterp, norm(len_annulus_min_initial - len_annulus_new_reinterp)); 
        end 
        
    end 
    
    valve_with_reference.leaflets(1).X = X;  


    
    debug = true; 
    if debug 
        valve_plot(valve_with_reference); 
    end 

end 



function Z = Out2(FUN,varargin)
    % Z = Out2(FUN,VARARGIN);
    %
    %	Provides the second output from the function
    [~,Z] = FUN(varargin{:});

end 