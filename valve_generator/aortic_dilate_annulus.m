function valve_with_reference = aortic_dilate_annulus(valve_with_reference)

    % dialation with preserved arc length 

    if isfield(valve_with_reference, 'dilation_dist')
        dilation_dist = valve_with_reference.dilation_dist; 
    else 
        error('must provide dilation_dist if valve.dilate_graft is true'); 
    end 

    leaflet = valve_with_reference.leaflets(1); 

    % compute annular length initial 
    [~, len_annulus_min_initial, len_annulus_each] = build_initial_fibers_aortic(leaflet, valve_with_reference); 
        
    % adjust radius 
    valve_with_reference.skeleton.r            = valve_with_reference.skeleton.r              + dilation_dist; 
    valve_with_reference.skeleton.r_of_z       = @(z) valve_with_reference.skeleton.r_of_z(z) + dilation_dist; 
    valve_with_reference.skeleton.r_commissure = valve_with_reference.skeleton.r_commissure   + dilation_dist; 
    % valve_with_reference.skeleton.r_co         = valve_with_reference.skeleton.r_co           + dilation_dist; 
    valve_with_reference.r                     = valve_with_reference.r                       + dilation_dist; 
    
    height_comm = valve_with_reference.skeleton.normal_height - valve_with_reference.skeleton.height_min_comm; 
    
    height_min_comm_override_initial_guess = valve_with_reference.skeleton.height_min_comm; 

    options = optimset('Display','off','TolFun',1e-16);

    % this has two return values 
    len_annulus_tmp = @(height_min_comm) build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm); 

    % hack to return only the next one 
    len_annulus = @(height_min_comm) Out2(len_annulus_tmp, height_min_comm); 
    
    len_annulus_minus_initial = @(height_min_comm) abs(len_annulus(height_min_comm) - len_annulus_min_initial); 

    % height_min_comm_new = fsolve(len_annulus_minus_initial,height_min_comm_override_initial_guess,options);     
    
    lower_bd = 0; 
    height_min_comm_new = fmincon(len_annulus_minus_initial,height_min_comm_override_initial_guess,[],[],[],[],lower_bd,[],[],options);     
    
    valve_with_reference.skeleton.height_min_comm = height_min_comm_new; 
    valve_with_reference.skeleton.normal_height = height_min_comm_new + height_comm; 
    
    [X, len_annulus_new] = build_initial_fibers_aortic(leaflet, valve_with_reference, height_min_comm_new); 

    tol = 1e-10;
    if norm(len_annulus_min_initial - len_annulus_new) > tol 
        warning('new annulus len not equal to previous after raising comm')
    end 
    
    
    redistribute_annulus_points = true; 
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
    
    % reset leaflet in data structure 
    valve_with_reference.leaflets(1).X = X;  
    
    
    

    % for j=1:j_max 
    %     for k=1:k_max 
    %        [th, r, z] = cart2pol(X(1,j,k), X(2,j,k), X(3,j,k));
    %        [x_tmp, y_tmp, z_tmp] = pol2cart(th,r + dilation_dist,z); 
    %        X(:,j,k) = [x_tmp; y_tmp; z_tmp]; 
    %     end 
    % end 
    % 
    % valve_with_reference.leaflets(1).X = X; 
    % 
    % valve_with_reference.skeleton.r = valve_with_reference.skeleton.r + dilation_dist; 
    % valve_with_reference.skeleton.r_of_z = @(z) valve_with_reference.skeleton.r_of_z(z) + dilation_dist; 

    debug = false; 
    if debug 
    %     figure; 
    % 
    %     x_component = squeeze(X(1,:,:)); 
    %     y_component = squeeze(X(2,:,:)); 
    %     z_component = squeeze(X(3,:,:)); 
    % 
    %     width = 1.0; 
    %     surf(x_component, y_component, z_component, 'LineWidth',width);
    %     axis equal 
    %     axis auto 

        valve_plot(valve_with_reference); 

    end 

end 


function Z = Out2(FUN,varargin)
    % Z = Out2(FUN,VARARGIN);
    %
    %	Provides the second output from the function
    [~,Z] = FUN(varargin{:});

end 