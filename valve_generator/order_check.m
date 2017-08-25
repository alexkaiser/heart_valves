function order_check(path, N_values)
    % 
    % 
    % 
    % 



    iterations = length(N_values) - 1; 

    diff_1    = zeros(iterations,1); 
    diff_2    = zeros(iterations,1); 
    diff_inf  = zeros(iterations,1); 

    diff_1_ring    = zeros(iterations,1); 
    diff_2_ring    = zeros(iterations,1); 
    diff_inf_ring  = zeros(iterations,1); 
    
    bilinear_interp = true; 

    if bilinear_interp 
        'bilinear interpolation (subject to interpolation error)'
    else 
        'injection interpolation (off by up to O(dx))'
    end 


    for i = 1:iterations 

        N_coarse = N_values(i); 
        data_coarse = sprintf('%s/mitral_tree_%d_final_data', path, N_coarse); 
        load(data_coarse) 

        valve_coarse       = valve; 
        leaflet_coarse     = valve_coarse.leaflets(1); 
        X_coarse           = leaflet_coarse.X; 
        j_max_anterior_coarse = max(leaflet_coarse.j_range_anterior); 
        j_max_coarse       = leaflet_coarse.j_max; 
        k_max_coarse       = leaflet_coarse.k_max; 
        du_coarse          = leaflet_coarse.du; 
        is_internal_coarse = leaflet_coarse.is_internal;

        N_fine = 2*N_coarse; 
        data_fine = sprintf('%s/mitral_tree_%d_final_data', path, N_fine); 
        load (data_fine) 


        valve_fine       = valve; 
        leaflet_fine     = valve_fine.leaflets(1); 
        X_fine           = leaflet_fine.X; 
        j_max_anterior_fine = max(leaflet_fine.j_range_anterior); 
        du_fine          = leaflet_fine.du; 
        j_max_fine       = leaflet_fine.j_max; 
        k_max_fine       = leaflet_fine.k_max; 
        is_internal_fine = leaflet_fine.is_internal; 
        is_bc_fine       = leaflet_fine.is_bc;    

        points_written_coarse = 0; 
        points_written_fine   = 0; 

        % points on X_coarse that have four valid neighbors 
        % on fine mesh for comparison 
        X_coarse_valid = zeros(size(X_coarse)); 

        % fine mesh interpolated to align with coarse mesh 
        X_fine_interp  = zeros(size(X_coarse)); 

        
        % valve ring point (for sanity checking)
        X_ring_coarse = zeros(size(X_coarse(:,:,k_max_coarse))); 
        X_ring_interp = zeros(size(X_ring_coarse)); 
        
        
        % for all internal points on coarse mesh 
        % find coordinates in u,v plane 
        % interpolate fine mesh to this location 
        % with bilinear interpolation 
        % only do so if all four points on the fine mesh are valid 
        for j=1:j_max_anterior_coarse 
        %for j=1:j_max_coarse
            for k=1:k_max_coarse
                
                % because of weird inclusive/exclusive angle placement 
                % anterior and posterior are every so slightly different here 
                
                if (j <= j_max_anterior_coarse)
                    N_coarse_anterior = j_max_anterior_coarse; 
                    
                    u = .5 * (j-1) / (N_coarse_anterior-1); 
                    v = .5 * (k-1) / (N_coarse_anterior-1); 
                else 
                    error('fix me!')
                    u = (j-1) / (N_coarse+1);
                    v = (k-1) / (N_coarse-1);                    
                end 
                
                N_fine_anterior = j_max_anterior_fine; 
                j_fine_below = floor(2 * u * (N_fine_anterior-1)) + 1; 
                k_fine_below = floor(2 * v * (N_fine_anterior-1)) + 1; 
                
                if is_internal_coarse(j,k)

                    % make sure that we are looking at a valid index 
                    if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && ... 
                       (0 < k_fine_below) && (k_fine_below <= k_max_fine) && ... 
                       (is_internal_fine(j_fine_below, k_fine_below) || is_bc_fine(j_fine_below, k_fine_below)) 

                        if bilinear_interp 

                            j_nbr_tmp = j_fine_below + 1; 
                            k_nbr_tmp = k_fine_below + 1; 

                            % check if have valid neighbor
                            % and get indices reduced for periodicity etc  
                            [valid j_fine_above k_fine_above] = get_indices(leaflet_fine, j_fine_below, k_fine_below, j_nbr_tmp, k_nbr_tmp);  

                            % above function handles periodicity and off to the side errors
                            % but does not tell if 
                            % because springs are owned by mimium index this is not checked correctly 
                            % and we need to make sure that all four locations are still valid 
                            if valid
                                if ~(is_internal_fine(j_fine_above, k_fine_below) || is_bc_fine(j_fine_above, k_fine_below)) 
                                    valid = false; 
                                end 
                                if ~(is_internal_fine(j_fine_below, k_fine_above) || is_bc_fine(j_fine_below, k_fine_above)) 
                                    valid = false; 
                                end
                                if ~(is_internal_fine(j_fine_above, k_fine_above) || is_bc_fine(j_fine_above, k_fine_above)) 
                                    valid = false; 
                                end
                            end 

                            if valid 

                                % have all necessary data to compare norms at this point 
                                % copy data into coarse 
                                X_coarse_valid(:,j,k) = X_coarse(:,j,k); 
                                points_written_coarse = points_written_coarse + 1; 

                                if (any(isnan(X_coarse_valid(:,j,k))) || any(X_coarse_valid(:,j,k) == [0;0;0]))
                                    error('nan or zero problems in coarse')
                                end 

                                % compute fractions of distances, which are the interpolation coefficients 
                                % this is whatever remains after taking floor 

                                % temp weight to only use bottom 
                                % this should still congeverge just with no particular order 


                                s = (u/du_fine) - (j_fine_below-1); 
                                t = (v/du_fine) - (k_fine_below-1); 

                                X_fine_interp(:,j,k) = (1-t) * ((1-s) * X_fine(:,j_fine_below,k_fine_below) + s*X_fine(:,j_fine_above,k_fine_below)) + ... 
                                                          t  * ((1-s) * X_fine(:,j_fine_below,k_fine_above) + s*X_fine(:,j_fine_above,k_fine_above)); 


                                if (any(isnan(X_fine_interp(:,j,k))) || any(X_fine_interp(:,j,k) == [0;0;0]))
                                    error('nan or zero problems in fine')
                                end 

                                points_written_fine = points_written_fine + 1;                          
                                
                            end 

                        else 

                            % injection, feels less rigorous in spirit 
                            % but really unclear if this is worse 

                            X_coarse_valid(:,j,k) = X_coarse(:,j,k); 
                            points_written_coarse = points_written_coarse + 1; 

                            if (any(isnan(X_coarse_valid(:,j,k))) || any(X_coarse_valid(1:2,j,k) == [0;0]))
                                % z component allowed to be zero
                                error('nan or zero problems in coarse')
                            end 

                            X_fine_interp(:,j,k) = X_fine(:,j_fine_below,k_fine_below); 

                            if (any(isnan(X_fine_interp(:,j,k))) || any(X_fine_interp(1:2,j,k) == [0;0]))
                                % z component allowed to be zero
                                error('nan or zero problems in fine')
                            end 

                            points_written_fine = points_written_fine + 1;                    

                        end 

                    end 
                end 
                
                if j==16 && (k == k_max_coarse)
                    'stop'
                end 
                
                % check at the ring 
                if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && (k == k_max_coarse)
                    s = (u/du_fine) - (j_fine_below-1); 
                    X_ring_coarse(:,j) = X_coarse(:,j,k); 
                    X_ring_interp(:,j) = (1-s) * X_fine(:,j_fine_below,k_max_fine) + s*X_fine(:,j_fine_above,k_max_fine); 
                end   
                
                
            end 
        end 

    %     % very simple sanity check 
    %     if ~all(find(X_coarse_valid(:)) == find(X_fine_interp(:)))
    %         warning('nonzero indices do not match, check results')
    %     end 

        diff_one_components = zeros(j_max_coarse,k_max_coarse); 
        diff_two_components = zeros(j_max_coarse,k_max_coarse); 
        diff_inf_components = zeros(j_max_coarse,k_max_coarse); 

        for j=1:j_max_coarse
            for k=1:k_max_coarse

                diff_one_components(j,k) = norm(X_coarse_valid(:,j,k) - X_fine_interp(:,j,k),1); 
                diff_two_components(j,k) = norm(X_coarse_valid(:,j,k) - X_fine_interp(:,j,k),2); 
                diff_inf_components(j,k) = norm(X_coarse_valid(:,j,k) - X_fine_interp(:,j,k),'inf'); 
            end 
        end 

        % these require a transpose to plot properly 
        diff_one_components = diff_one_components'; 
        diff_two_components = diff_two_components'; 
        diff_inf_components = diff_inf_components'; 

        u_mesh = (1:j_max_coarse) / N_coarse; 
        v_mesh = (1:k_max_coarse) / N_coarse; 

        
        
        % interpolation method check on ring 
        diff_one_ring_components = zeros(j_max_coarse,1); 
        diff_two_ring_components = zeros(j_max_coarse,1); 
        diff_inf_ring_components = zeros(j_max_coarse,1); 

        for j=1:j_max_coarse
            diff_one_ring_components(j) = norm(X_ring_coarse(:,j) - X_ring_interp(:,j),1); 
            diff_two_ring_components(j) = norm(X_ring_coarse(:,j) - X_ring_interp(:,j),2); 
            diff_inf_ring_components(j) = norm(X_ring_coarse(:,j) - X_ring_interp(:,j),'inf'); 
        end  
        
        differences = X_ring_coarse(:) - X_ring_interp(:);  

        % du measure for single integral 
        diff_1_ring(i)    = du_coarse       * norm(differences, 1); 

        % on 2 norm gets a root on the measure 
        diff_2_ring(i)    = sqrt(du_coarse) * norm(differences, 2); 

        % no weight on sup norm version 
        diff_inf_ring(i)  =                   norm(differences, 'inf'); 
        
        
        
        
    %     figure; 
    %     diff_one_components(diff_one_components == 0) = nan;   
    %     surf(diff_one_components); 
    %     axis tight; 
    %     figure; 
    %     diff_two_components(diff_two_components == 0) = nan; 
    %     surf(diff_two_components); 
    %     axis tight; 

        x_min = 0; 
        x_max = 1; 
        y_min = .2;
        y_max = max(v_mesh);
        z_min = 0;
        z_max = max(max(diff_inf_components)); 

        figure; 
        diff_inf_components(diff_inf_components == 0) = nan; 
        surf(u_mesh, v_mesh, diff_inf_components); 
        axis([x_min x_max y_min y_max z_min z_max]); 
        xlabel('u')
        ylabel('v')


        figure; 
        diff_inf_components(diff_inf_components == 0) = nan; 
        surf(u_mesh, v_mesh, diff_inf_components,'EdgeColor','None'); 
        axis([x_min x_max y_min y_max z_min z_max]); 
        view(2);  
        xlabel('u')
        ylabel('v')

%         figure 
%         pcolor(u_mesh, v_mesh, diff_inf_components)
%         shading flat 
%         axis([x_min x_max y_min y_max]);
%         xlabel('u')
%         ylabel('v')



        % For results to be comparable across vectors of different sizes
        % Need to scale norms so that they are double integrals 

        differences = X_coarse_valid(:) - X_fine_interp(:);  

        % du^2 measure for double integral 
        diff_1(i)    = du_coarse*du_coarse * norm(differences, 1); 

        % on 2 norm gets a root on the measure 
        diff_2(i)    = du_coarse           * norm(differences, 2); 

        % no weight on sup norm version 
        diff_inf(i)  =                       norm(differences, 'inf'); 

    end
    
    diff_1_ring 
    diff_2_ring
    diff_inf_ring 


    fprintf('Comparisons of N,2N and 2N,4N points\n'); 
    N = N_values(1); 
    for i = 1:iterations-1

        order_1    = diff_1_ring(i) / diff_1_ring(i+1); 
        order_2    = diff_2_ring(i) / diff_2_ring(i+1);  
        order_inf  = diff_inf_ring(i) / diff_inf_ring(i+1);  

        fprintf('%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N,  order_1, order_2, order_inf);

        N = 2*N;
    end
    
    
    
    
    diff_1
    diff_2 
    diff_inf 


    fprintf('Comparisons of N,2N and 2N,4N points\n'); 
    N = N_values(1); 
    for i = 1:iterations-1

        order_1    = diff_1(i) / diff_1(i+1); 
        order_2    = diff_2(i) / diff_2(i+1);  
        order_inf  = diff_inf(i) / diff_inf(i+1);  

        fprintf('%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N,  order_1, order_2, order_inf);

        N = 2*N;
    end

end 



