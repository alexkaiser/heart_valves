function order_check(path, N_values, suffix_name)
    % 
    % 
    % 
    % 

    ring_convergence_sanity_check = false; 
    front_view_plots = true; 

    iterations = length(N_values) - 1; 

    diff_1    = zeros(iterations,1); 
    diff_2    = zeros(iterations,1); 
    diff_inf  = zeros(iterations,1); 

    diff_1_anterior     = zeros(iterations,1); 
    diff_2_anterior     = zeros(iterations,1); 
    diff_inf_anterior   = zeros(iterations,1);
    
    diff_1_posterior    = zeros(iterations,1); 
    diff_2_posterior    = zeros(iterations,1); 
    diff_inf_posterior  = zeros(iterations,1);
    
    diff_1_ring    = zeros(iterations,1); 
    diff_2_ring    = zeros(iterations,1); 
    diff_inf_ring  = zeros(iterations,1); 
    
    newton_error_difference_equations = zeros(length(N_values),1); 
    
    bilinear_interp = true; 

    if bilinear_interp 
        fprintf('bilinear interpolation (subject to interpolation error)\n'); 
    else 
        fprintf('injection interpolation (off by up to O(dx))\n'); 
    end 


    for i = 1:iterations 

        N_coarse = N_values(i)
        data_coarse = sprintf('%s/mitral_tree_%d_final_data', path, N_coarse); 
        load(data_coarse, 'valve'); 

        valve_coarse       = valve; 
        leaflet_coarse     = valve_coarse.leaflets(1); 
        X_coarse           = leaflet_coarse.X; 
        j_max_anterior_coarse = max(leaflet_coarse.j_range_anterior); 
        j_max_coarse       = leaflet_coarse.j_max; 
        k_max_coarse       = leaflet_coarse.k_max; 
        du_coarse          = leaflet_coarse.du; 
        is_internal_coarse = leaflet_coarse.is_internal; 

        newton_error_difference_equations(i) = total_global_err(leaflet_coarse); 
        
        if front_view_plots 
            
            fig = figure; 
            leaflet_coarse_anterior_only = leaflet_coarse; 
            
            % nan mask posterior 
            leaflet_coarse_anterior_only.X(:,leaflet_coarse.j_range_posterior,:) = nan; 
            
            % junk hack, say there are only two trees 
            leaflet_coarse_anterior_only.num_trees = 2; 
            
            surf_plot(leaflet_coarse_anterior_only, fig); 
            axis equal; 
            printfig(fig, sprintf('diag_view_anterior_%d_%s.eps', N_coarse, suffix_name));
            
            fig = figure; 
            leaflet_coarse_posterior_only = leaflet_coarse; 
            
            % nan mask posterior 
            leaflet_coarse_posterior_only.X(:,leaflet_coarse.j_range_anterior,:) = nan; 
            
            % junk hack, say there are only two trees 
            leaflet_coarse_posterior_only.num_trees = 6; 
            leaflet_coarse_posterior_only.chordae = leaflet_coarse_posterior_only.chordae(3:8); 
            
            surf_plot(leaflet_coarse_posterior_only, fig); 
            axis equal; 
            printfig(fig, sprintf('diag_view_posterior_%d_%s.eps',N_coarse, suffix_name));
            
        end 
        
        if (~isempty(leaflet_coarse.j_range_left_comm)) || ...
           (~isempty(leaflet_coarse.j_range_right_comm))
            error('Order check not implemented with commissural leaflets'); 
        end 
        
%         if i == 1
%             fig = figure; 
%             valve_plot(valve_coarse, fig);
%         end 
        
        N_fine = 2*N_coarse; 
        data_fine = sprintf('%s/mitral_tree_%d_final_data', path, N_fine); 
        load (data_fine, 'valve') 


        valve_fine       = valve; 
        leaflet_fine     = valve_fine.leaflets(1); 
        X_fine           = leaflet_fine.X; 
        j_max_anterior_fine = max(leaflet_fine.j_range_anterior); 
        du_fine          = leaflet_fine.du; 
        j_max_fine       = leaflet_fine.j_max; 
        k_max_fine       = leaflet_fine.k_max; 
        is_internal_fine = leaflet_fine.is_internal; 
        is_bc_fine       = leaflet_fine.is_bc;    

        % final iteration 
        if (i==iterations)
            N_fine
            newton_error_difference_equations(i+1) = total_global_err(leaflet_fine); 
            
            if front_view_plots 
                fig = figure; 
                leaflet_fine_anterior_only = leaflet_fine; 

                % nan mask posterior 
                leaflet_fine_anterior_only.X(:,leaflet_fine.j_range_posterior,:) = nan; 

                % junk hack, say there are only two trees 
                leaflet_fine_anterior_only.num_trees = 2; 

                surf_plot(leaflet_fine_anterior_only, fig); 
                axis equal; 
                printfig(fig, sprintf('diag_view__anterior_%d_%s.eps',N_fine, suffix_name));
                
                fig = figure; 
                leaflet_fine_posterior_only = leaflet_fine; 

                % nan mask posterior 
                leaflet_fine_posterior_only.X(:,leaflet_fine.j_range_anterior,:) = nan; 

                % junk hack, say there are only two trees 
                leaflet_fine_posterior_only.num_trees = 6; 
                leaflet_fine_posterior_only.chordae = leaflet_fine_posterior_only.chordae(3:8); 
                
                surf_plot(leaflet_fine_posterior_only, fig); 
                axis equal; 
                printfig(fig, sprintf('diag_view_posterior_%d_%s.eps',N_fine, suffix_name));
                
            end 
            
        end 
        
%         fig = figure; 
%         valve_plot(valve_fine, fig);
        
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
        
        
        u_mesh = zeros(j_max_coarse,1); 
        v_mesh = zeros(k_max_coarse,1); 
        
        % for all internal points on coarse mesh 
        % find coordinates in u,v plane 
        % interpolate fine mesh to this location 
        % with bilinear interpolation 
        % only do so if all four points on the fine mesh are valid 
        %for j=1:j_max_anterior_coarse 
        for j=1:j_max_coarse
        %for j=(j_max_anterior_coarse+1):j_max_coarse 
            for k=1:k_max_coarse
                
                % because of weird inclusive/exclusive angle placement 
                % anterior and posterior are every so slightly different here 
                N_coarse_anterior = j_max_anterior_coarse; 
                N_fine_anterior = j_max_anterior_fine; 
                
                N_coarse_posterior = N_coarse - N_coarse_anterior; 
                N_fine_posterior = N_fine - N_fine_anterior;  
                
                if (j <= j_max_anterior_coarse)
                    
                    u = .5 * (j-1) / (N_coarse_anterior-1);
                    j_fine_below = floor(2 * u * (N_fine_anterior-1)) + 1;
                    
                    % compute fractions of distances, which are the interpolation coefficients 
                    s = 2*u*(N_fine_anterior-1) - (j_fine_below - 1); 
                    
                else
                    
                    % index in the posterior range 
                    j_reduced = j - j_max_anterior_coarse;
                    
                    % posterior is exclusive of ends, take plus one on N
                    u = .5 + .5 * j_reduced / (N_coarse_posterior + 1);  
                    
                    j_fine_below_reduced = floor(2*(u-.5) * (N_fine_posterior+1)); 
                    
                    
                    
                    % compute fractions of distances, which are the interpolation coefficients 
                    s = (2*(u-.5) * (N_fine_posterior + 1)) - (j_fine_below_reduced); 
                    
                    j_fine_below = j_fine_below_reduced + j_max_anterior_fine; 
                    
                    u_fine_check = .5 + .5 * j_fine_below_reduced / (N_fine_posterior + 1);  
                end 
                
                if ~bilinear_interp
                    s = 0; 
                end 
                                
                % v is independent of anterior or posterior 
                v = 1 - (k_max_coarse - k) * du_coarse; 
                                
                % new fancy... 
                k_fine_below = (v - 1)/du_fine + k_max_fine; 
                
                v_fine_check = 1 - (k_max_fine - k_fine_below) * du_fine; 
                
                % fprintf('k = %d, k_fine_below = %d, v = %f, v_fine_check = %f\n', k, k_fine_below, v, v_fine_check); 
                
                
                % for plotting 
                u_mesh(j) = u; 
                v_mesh(k) = v; 
                
                
                if k_fine_below ~= floor(k_fine_below)
                    warning('something off on k indexing, not rounding correctly')
                end 
                
                if abs(v - v_fine_check) > eps 
                    error('v and v_fine_check not equal to tolerance')
                end 
                
                if is_internal_coarse(j,k)

                    % make sure that we are looking at a valid index 
                    if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && ... 
                       (0 < k_fine_below) && (k_fine_below <= k_max_fine) && ... 
                       (is_internal_fine(j_fine_below, k_fine_below) || is_bc_fine(j_fine_below, k_fine_below)) 

                        if bilinear_interp 

                            j_nbr_tmp = j_fine_below + 1; 
                            k_nbr_tmp = k_fine_below; 

                            % check if have valid neighbor
                            % and get indices reduced for periodicity etc  
                            [valid j_fine_above] = get_indices(leaflet_fine, j_fine_below, k_fine_below, j_nbr_tmp, k_nbr_tmp);  

                            % above function handles periodicity and off to the side errors
                            % but does not tell if 
                            % because springs are owned by mimium index this is not checked correctly 
                            % and we need to make sure that all four locations are still valid 
                            if valid
                                if ~(is_internal_fine(j_fine_above, k_fine_below) || is_bc_fine(j_fine_above, k_fine_below)) 
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

                                if s == 0 
                                    % edge case, s==0 may be the last point (and there is no point aboce to be read)
                                    X_fine_interp(:,j,k) = X_fine(:,j_fine_below,k_fine_below);
                                else 
                                    X_fine_interp(:,j,k) = (1-s) * X_fine(:,j_fine_below,k_fine_below) + s*X_fine(:,j_fine_above,k_fine_below);
                                end 
                                
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
                    'stop'; 
                end 
                
                % check at the ring to make sure checker gets desired order 
                % piecewise linear interpolant (with exact application of boundary conditions)
                % should give secon order on ring 
                if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && (k == k_max_coarse) 
                    
                    % fprintf('j = %d, j_fine_below = %d, u = %f, s = %f\n', j, j_fine_below, u, s); 
                    X_ring_coarse(:,j) = X_coarse(:,j,k); 

                    if s == 0
                        % edge case, s==0 may be the last point (and there is no point above to be read)
                        X_ring_interp(:,j) = X_fine(:,j_fine_below,k_max_fine); 
                    else 
                        X_ring_interp(:,j) = (1-s) * X_fine(:,j_fine_below,k_max_fine) + s*X_fine(:,j_fine_above,k_max_fine); 
                    end 
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
        y_min = .7; 
        y_max = max(v_mesh);
        z_min = 0;
        z_max = max(max(diff_two_components)) + .01; 

%         figure; 
%         diff_inf_components(diff_inf_components == 0) = nan; 
%         surf(u_mesh, v_mesh, diff_inf_components); 
%         axis([x_min x_max y_min y_max z_min z_max]); 
%         xlabel('u')
%         ylabel('v')

        if N_fine == 1024

            fig = figure; 
            set(gcf,'Renderer','Zbuffer')
            diff_two_components(diff_two_components == 0) = nan; 
            surf(u_mesh, v_mesh, diff_two_components,'EdgeColor','None'); 
            set(fig, 'Position', [100, 100, 1000, 500])
            set(fig,'PaperPositionMode','auto')
            axis equal
            axis([x_min x_max y_min y_max z_min z_max]);
            view(2);  
            xlabel('u')
            ylabel('v')
            cb = colorbar('SouthOutside'); 
            set(cb,'position',[0.13   0.16   0.35   0.02])
            
            cb.Label.String = 'cm'; 
            %set(cb,'position',[0.13   0.22   0.78   0.02])

            
            printfig(fig,sprintf('static_convergence_psuedocolor_difference_%s', suffix_name)); 
            
            
            fig = figure; 
            diff_two_components(diff_two_components == 0) = nan; 
            surf(u_mesh, v_mesh, diff_two_components, 'EdgeColor','None');  
            set(fig, 'Position', [100, 100, 1000, 500])
            set(fig,'PaperPositionMode','auto')
            axis equal
            axis([x_min x_max y_min y_max z_min z_max]);
            grid off 
            axis off 
            view(-53,14);  
            %cb = colorbar('location', 'EastOutside'); 
            %set(cb,'position',[0.8    0.38    0.0446    0.3])
            %xlabel('u')
            %ylabel('v')
            printfig(fig, sprintf('static_convergence_surface_difference_%s', suffix_name)); 
        
        end 

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

        % anterior only for diagnosting problems 
        differences_anterior = X_coarse_valid(:,leaflet_coarse.j_range_anterior,:) - X_fine_interp(:,leaflet_coarse.j_range_anterior,:);  
        differences_anterior = differences_anterior(:); 

        % du^2 measure for double integral 
        diff_1_anterior(i)    = du_coarse*du_coarse * norm(differences_anterior, 1); 

        % on 2 norm gets a root on the measure 
        diff_2_anterior(i)    = du_coarse           * norm(differences_anterior, 2); 

        % no weight on sup norm version 
        diff_inf_anterior(i)  =                       norm(differences_anterior, 'inf'); 
        
        % posterior only for diagnosting problems 
        differences_posterior = X_coarse_valid(:,leaflet_coarse.j_range_posterior,:) - X_fine_interp(:,leaflet_coarse.j_range_posterior,:);  
        differences_posterior = differences_posterior(:); 

        % du^2 measure for double integral 
        diff_1_posterior(i)    = du_coarse*du_coarse * norm(differences_posterior, 1); 

        % on 2 norm gets a root on the measure 
        diff_2_posterior(i)    = du_coarse           * norm(differences_posterior, 2); 

        % no weight on sup norm version 
        diff_inf_posterior(i)  =                       norm(differences_posterior, 'inf'); 
    end
    
    if ring_convergence_sanity_check
        diff_1_ring 
        diff_2_ring
        diff_inf_ring 


        fprintf('Comparisons of N,2N and 2N,4N points\n'); 
        N = N_values(1); 
        for i = 1:iterations-1

            order_1    = diff_1_ring(i) / diff_1_ring(i+1); 
            order_2    = diff_2_ring(i) / diff_2_ring(i+1);  
            order_inf  = diff_inf_ring(i) / diff_inf_ring(i+1);  

            fprintf('%d,%d,%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N, 2*N, 4*N, order_1, order_2, order_inf);

            N = 2*N;
        end
    end 
    

    fprintf('Difference equations 2 norm (i.e. Newton solve two norm error):\n')
    N = N_values(1); 
    for i=1:length(N_values)
        fprintf('%d &  %e \\\\  \n \\hline \n ', N, newton_error_difference_equations(i)); 
        N = 2*N;
    end 
    
    
    
    fprintf('Comparisons of N,2N points\n'); 
    N = N_values(1); 
    for i = 1:iterations

        fprintf('%d,%d & %.2e & %.2e  &  %.2e    \\\\  \n \\hline \n ',  N, 2*N,   diff_1(i), diff_2(i), diff_inf(i));

        N = 2*N;
    end
    
    
    fprintf('Comparisons of N,2N and 2N,4N points\n'); 
    N = N_values(1); 
    for i = 1:iterations-1

        order_1    = diff_1(i) / diff_1(i+1); 
        order_2    = diff_2(i) / diff_2(i+1);  
        order_inf  = diff_inf(i) / diff_inf(i+1);  

        fprintf('%d,%d,%d & %.2f & %.2f  &  %.2f    \\\\  \n \\hline \n ',  N, 2*N, 4*N,   order_1, order_2, order_inf);

        N = 2*N;
    end

    
    
    fprintf('Anterior, Comparisons of N,2N points\n'); 
    N = N_values(1); 
    for i = 1:iterations

        fprintf('%d,%d & %.2e & %.2e  &  %.2e    \\\\  \n \\hline \n ',  N, 2*N,   diff_1_anterior(i), diff_2_anterior(i), diff_inf_anterior(i));

        N = 2*N;
    end
    
    
    fprintf('Anterior, Comparisons of N,2N and 2N,4N points\n'); 
    N = N_values(1); 
    for i = 1:iterations-1

        order_1    = diff_1_anterior(i) / diff_1_anterior(i+1); 
        order_2    = diff_2_anterior(i) / diff_2_anterior(i+1);  
        order_inf  = diff_inf_anterior(i) / diff_inf_anterior(i+1);  

        fprintf('%d,%d,%d & %.2f & %.2f  &  %.2f    \\\\  \n \\hline \n ',  N, 2*N, 4*N,   order_1, order_2, order_inf);

        N = 2*N;
    end

    
    
    fprintf('posterior, Comparisons of N,2N points\n'); 
    N = N_values(1); 
    for i = 1:iterations

        fprintf('%d,%d & %.2e & %.2e  &  %.2e    \\\\  \n \\hline \n ',  N, 2*N,   diff_1_posterior(i), diff_2_posterior(i), diff_inf_posterior(i));

        N = 2*N;
    end
    
    
    fprintf('Posterior, Comparisons of N,2N and 2N,4N points\n'); 
    N = N_values(1); 
    for i = 1:iterations-1

        order_1    = diff_1_posterior(i) / diff_1_posterior(i+1); 
        order_2    = diff_2_posterior(i) / diff_2_posterior(i+1);  
        order_inf  = diff_inf_posterior(i) / diff_inf_posterior(i+1);  

        fprintf('%d,%d,%d & %.2f & %.2f  &  %.2f    \\\\  \n \\hline \n ',  N, 2*N, 4*N,   order_1, order_2, order_inf);

        N = 2*N;
    end
    
end 



