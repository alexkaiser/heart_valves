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

%     if N_leaflets ~= 2
%         error("N_leaflets must be equal to 2")
%     end 

    k_range = 2:k_max; 

    if ~exist('extra_stretch_radial', 'var')
        extra_stretch_radial = 1.0; 
    end

    if ~exist('extra_stretch_circ', 'var')
        extra_stretch_circ = 1.0; 
    end

    % search down in power to find as flat a leaflet as possible by the comms 
    % but without intersection 
    % power_search_list = [4,3,2,1.75,1.5,1.25,1];
    power_search_list = [2,1.75,1.5,1.25,1];
    % power_search_list = [4];
    
    pass_y_gt_0_constraint = false; 
    
    for power = power_search_list
    
        X = leaflet.X; 
        
        % tol change 
        tol_rms_err_strain = 1e-10; 

        % start with initial free edge 
        % commissure points stay fixed 
        % free_edge_interp_points = X(:,:,k_max); 

        % give it plenty 
        n_iterations = 2000; 

        debug_lengths = false; 

        % for debug info 
        [free_edge_length_single_loaded, free_edge_length_single_rest] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths); 

        % use the intercomm radius
        r = leaflet.skeleton.r_commissure; 
        % find radius of annulus
        % r = norm(leaflet.X(1:2,1,1));  


        % find desired coefficient for initial curve before iteration 
        % this is the parameter to search over 
        % y_max_from_center = 0.9; 

        free_edge_len_minus_rest = @(y_max) abs(run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max, power) - free_edge_length_single_rest * extra_stretch_circ); 
        

        options = optimset('Display','off','TolFun',1e-16);

        % compute original midpoint of chord
        % in xy plane 
        % bounds relative to this position 
        chord_midpoint = 0.5 * (X(1:2,1,k_max) + X(1:2,j_max,k_max));
        norm_midpoint = norm(chord_midpoint);
        

        % do not go lower than threshold 
        % can be close in the middle 
        % originally tuned for r = 1.25 cm valve, scale relative to this 
        r_basic = 1.25; 
                
        if N_leaflets == 2
            % for two leaflets away from center more 
            y_max_from_center_min = (r/r_basic) * 0.4;
            % and keep away from from the outside 
            y_max_from_center_max_thresh = (r/r_basic) * (r_basic - 0.4);
            % place on radius in two leaflet case 
            y_max_from_center_initial_guess = 1.0 * (r/r_basic); 
        else 
            % 3 or more leaflets 
            
            % allow to go fairly far in 
            y_max_from_center_min = (r/r_basic) * 0.2;
            % and keep away from from the outside 
            y_max_from_center_max_thresh = (r/r_basic) * (r_basic - 0.4);
            
            % start at center 
            y_max_from_center_initial_guess = -r;
        end 

        min_bound = y_max_from_center_min - norm_midpoint;
        max_bound = y_max_from_center_max_thresh - norm_midpoint;
        
        use_fsolve = false; 
        use_fmincon = true; 

        if use_fsolve 
            y_max_from_center = fsolve(free_edge_len_minus_rest,y_max_from_center_initial_guess,options); 
            y_max_from_center = max(y_max_from_center, y_max_from_center_min); 
        elseif use_fmincon
%             y_max_from_center = fmincon(free_edge_len_minus_rest, y_max_from_center_initial_guess,[],[],[],[],y_max_from_center_min,y_max_from_center_max_thresh,[],options)
            y_max_from_center = fmincon(free_edge_len_minus_rest, y_max_from_center_initial_guess,[],[],[],[],min_bound,max_bound,[],options)
        else 
            y_max_from_center = r/2; 
        end 

        r
        y_max_from_center

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


        % compute the resulting values 
        [free_edge_len_temp, free_edge_interp_points] = run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max_from_center, power);

        pass_equal_strain = false; 

        for it = 1:n_iterations

            for j=2:N_each
                for k=k_range

                    ring_point = X(:,j,1); 

                    % total radial rest length of this radial fiber 
                    total_rest_length = sum(R_v(j, 1:(k-1))); 

                    tangent = (free_edge_interp_points(:,j) - ring_point); 
                    tangent = tangent / norm(tangent); 

                    % based on the rest length 
                    X(:,j,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

                end 
            end 



            % get free edge lengths 
            % 'after sin^2 interpolation'
            [free_edge_length_single_loaded, free_edge_length_single_rest, portion_of_current_edge, portion_of_free_edge] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
            % free_edge_length_single_loaded, free_edge_length_single_rest; 

            X_free_edge_leaflet_1 = X(:,:,k_max); 

            % spacing of points as fraction of arc length 
            interp_idx_free_edge = [0; portion_of_current_edge];  

            % interpolate as fraction of rest length 
            free_edge_interp_points_respaced = interp1(interp_idx_free_edge, X_free_edge_leaflet_1', portion_of_free_edge(1:N_each-1))'; 

            free_edge_interp_points(:, 2:N_each) = free_edge_interp_points_respaced; 

            X(:,:,k_max) = free_edge_interp_points; 

            % 'after respacing free edge interp points'
            [free_edge_length_single_loaded, free_edge_length_single_rest, ~, ~, strains] = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u, debug_lengths);
            % free_edge_length_single_loaded, free_edge_length_single_rest

            mean_strain = mean(strains); 
            rms_err_strain = sqrt(sum((strains - mean_strain).^2)); 

            if rms_err_strain < tol_rms_err_strain
                % fprintf('Exiting on it %d with mean strain %f and rms error on strain %e\n', it, mean_strain, rms_err_strain); 
                pass_equal_strain = true; 
                break
            end 

        end % strain equalize loop 

        fprintf('Power %f through strain loop, pass_equal_strain = %d\n', power, pass_equal_strain);
        
        if N_leaflets == 2
            if all(free_edge_interp_points(2,2:N_each) > 0)
                pass_y_gt_0_constraint = true; 
                fprintf('Power %f passed nonzero y constraint\n', power);
                break; 
            else 
                % pass this constraint already false 
                % pass_y_gt_0_constraint = false; 
                fprintf('Found negative in free_edge_interp_points with power %f\n', power);
            end 
        end 
    
        if N_leaflets > 2
            if pass_equal_strain
                break; 
            end 
        end 
        
    end % power outer loop
    
   
    % cleanup 
    if ~pass_equal_strain 
        warning('Exiting on max it %d with mean strain %f and rms error on strain %e, not fully converged.\n', it, mean_strain, rms_err_strain')
    end
    
    if ~pass_y_gt_0_constraint
        warning('Exiting with power %f and negative values for y along free edge\n');
    end 
        
end 



