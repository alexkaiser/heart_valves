function free_edge_interp_points = find_free_edge_interp_points_true_bicuspid(leaflet, extra_stretch_radial, extra_stretch_circ)

    j_max  = leaflet.j_max; 
    k_max  = leaflet.k_max; 
    N_each = leaflet.N_each; 

    R_u    = leaflet.R_u; 
    R_v    = leaflet.R_v; 

    if isfield(leaflet, 'N_leaflets')
        N_leaflets = leaflet.N_leaflets; 
    else 
        error("must provide N_leaflets")
    end 

    if N_leaflets ~= 2
        error("N_leaflets must be equal to 2")
    end 

    k_range = 2:k_max; 

    if ~exist('extra_stretch_radial', 'var')
        extra_stretch_radial = 1.0; 
    end

    if ~exist('extra_stretch_circ', 'var')
        extra_stretch_circ = 1.0; 
    end
    
    X = leaflet.X; 



    % tol change 
    tol_rms_err_strain = 1e-10; 
    
    % start with initial free edge 
    % commissure points stay fixed 
    free_edge_interp_points = X(:,:,k_max); 
    
    % give it plenty 
    n_iterations = 2000; 
    
    debug_lengths = false; 
    
    % for debug info 
    [free_edge_length_single_loaded, free_edge_length_single_rest] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths); 
    
    % find radius of annulus
    r = norm(leaflet.X(1:2,1,1));  
    
    
    % find desired coefficient for initial curve before iteration 
    % this is the parameter to search over 
    % y_max_from_center = 0.9; 
    
    free_edge_len_minus_rest = @(y_max) abs(run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max) - free_edge_length_single_rest * extra_stretch_circ); 
    % free_edge_len_minus_rest = @(y_max) abs(run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max)); 
    

    options = optimset('Display','off','TolFun',1e-16);

    % 
    
    % do not go lower than threshold 
    % can be close in the middle 
    % originally tuned for r = 1.25 cm valve, scale relative to this 
    r_basic = 1.25; 
    y_max_from_center_min = (r/r_basic) * 0.1;
    % but keep farther from the outside 
    y_max_from_center_max_thresh = (r/r_basic) * (r_basic - 0.4);
    
    y_max_from_center_initial_guess = 1.0 * (r/r_basic); 
    
    use_fsolve = false; 
    use_fmincon = true; 
    
    if use_fsolve 
        y_max_from_center = fsolve(free_edge_len_minus_rest,y_max_from_center_initial_guess,options); 
        y_max_from_center = max(y_max_from_center, y_max_from_center_min); 
    elseif use_fmincon
        y_max_from_center = fmincon(free_edge_len_minus_rest, y_max_from_center_initial_guess,[],[],[],[],y_max_from_center_min,y_max_from_center_max_thresh,[],options);           
    else 
        y_max_from_center = r/2; 
    end 
    
%     y_max_from_center_fsolve = fsolve(free_edge_len_minus_rest,y_max_from_center_initial_guess,options)
%     y_max_from_center_fmincon_unconstrained = fmincon(free_edge_len_minus_rest, y_max_from_center_initial_guess,[],[],[],[],[],[],[],options) 
%     y_max_from_center_fmincon = fmincon(free_edge_len_minus_rest, y_max_from_center_initial_guess,[],[],[],[],y_max_from_center_min,y_max_from_center_max_thresh,[],options) 
%     
    
    debug_plots = false; 
    if debug_plots
        figure; 
        y_range = linspace(0,r,1000); 
        
        vals = zeros(size(y_range)); 
        
        for j=1:length(y_range)
            vals(j) = free_edge_len_minus_rest(y_range(j)); 
        end 
        
        plot(y_range, vals); 
        % xlabel('y', 'free_edge_len_minus_rest'); 
    end 
    

    
    % set initial version of the new leaflet position 
    % with free edge at given functional form (sin^2) 
    for comm_idx = 1:N_leaflets

        % point one internal of commissure to point that m
        % N_each is a power of two 
        min_idx = (comm_idx-1)*N_each;         

        for j=1:(N_each-1)

            ring_point = X(:,j + min_idx ,1); 

            th = atan2(ring_point(2), ring_point(1));  

            % make a litle 
            strained_len_total = extra_stretch_radial * sum(R_v(j + min_idx, :)); 

            % cm apart at middle 
            y_free_edge_end = y_max_from_center * sign(sin(th)) * sin(th)^2; 
            % y_free_edge_end = 0; 
            % this would put the two free edges exactly coinciding 

            % if using exact x 
            % then (y_diff^2 + height^2) = strained_len_total^2 
            % so height is given as 
            interp_height = sqrt(strained_len_total^2 - (ring_point(2) - y_free_edge_end)^2) ; 

            % comm_interp_point = [ring_point(1) ; (r/2) * sin(th); comm_prev(3)];                 
            free_edge_interp_points(:,j + min_idx) = [ring_point(1) ; y_free_edge_end; interp_height + ring_point(3)]; 

        end 
    end 

    
    pass = false; 
    
    for it = 1:n_iterations
                
        % sets the new leaflet position 
        % with free edge at given functional form (sin^2) 
        for comm_idx = 1:N_leaflets

            % point one internal of commissure to point that m
            % N_each is a power of two 
            min_idx = (comm_idx-1)*N_each;         

            for j=1:(N_each-1)
                for k=k_range

                    ring_point = X(:,j + min_idx ,1); 

                    % total radial rest length of this radial fiber 
                    total_rest_length = sum(R_v(j + min_idx, 1:(k-1))); 

                    tangent = (free_edge_interp_points(:,j) - ring_point); 
                    tangent = tangent / norm(tangent); 

                    % based on the rest length 
                    X(:,j + min_idx ,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

                end 
            end 

        end 
        
        % get free edge lengths 
        % 'after sin^2 interpolation'
        [free_edge_length_single_loaded, free_edge_length_single_rest, portion_of_current_edge, portion_of_free_edge] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
        % free_edge_length_single_loaded, free_edge_length_single_rest; 

        X_free_edge_leaflet_1_with_wrap = [X(:,j_max,k_max), X(:, 1:N_each,k_max)]; 

        % spacing of points as fraction of arc length 
        interp_idx_free_edge = [0; portion_of_current_edge];  

        % interpolate as fraction of rest length 
        free_edge_interp_points_respaced = interp1(interp_idx_free_edge, X_free_edge_leaflet_1_with_wrap', portion_of_free_edge(1:N_each-1))'; 

        free_edge_interp_points(:, (1:N_each-1)) = free_edge_interp_points_respaced; 

        % other leaflet by rotation 
        free_edge_interp_points(:, (N_each+1:j_max-1)) = rotation_matrix_z(pi) * free_edge_interp_points_respaced; 

        X(:,:,k_max) = free_edge_interp_points; 

        % 'after respacing free edge interp points'
        [free_edge_length_single_loaded, free_edge_length_single_rest, ~, ~, strains] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
        % free_edge_length_single_loaded, free_edge_length_single_rest
        
        
        mean_strain = mean(strains); 
        rms_err_strain = sqrt(sum((strains - mean_strain).^2)); 
        
        if rms_err_strain < tol_rms_err_strain
            % fprintf('Exiting on it %d with mean strain %f and rms error on strain %e\n', it, mean_strain, rms_err_strain); 
            pass = true; 
            break
        end 
        
    end 
    
    if ~pass 
        warning('Exiting on max it %d with mean strain %f and rms error on strain %e, not fully converged.\n', it, mean_strain, rms_err_strain')
    end 
    
    
end 


