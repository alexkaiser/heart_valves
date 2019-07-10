function order_check_aortic(path, N_values, suffix_name)
    % 
    % 
    % 
    % 

    % Copyright (c) 2019, Alexander D. Kaiser
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
    iterations = length(N_values) - 1; 

    diff_1    = zeros(iterations,1); 
    diff_2    = zeros(iterations,1); 
    diff_inf  = zeros(iterations,1); 

    newton_error_difference_equations = zeros(length(N_values),1); 

    for i = 1:iterations 

        N_coarse = N_values(i)
        data_coarse = sprintf('%s/aortic_%d_final_data', path, N_coarse); 
        load(data_coarse, 'valve'); 

        valve_coarse       = valve; 
        leaflet_coarse     = valve_coarse.leaflets(1); 
        X_coarse           = leaflet_coarse.X; 
        j_max_coarse       = leaflet_coarse.j_max; 
        k_max_coarse       = leaflet_coarse.k_max; 
        du_coarse          = leaflet_coarse.du; 

        newton_error_difference_equations(i) = total_global_err(leaflet_coarse); 
        
        N_fine = 2*N_coarse; 
        data_fine = sprintf('%s/aortic_%d_final_data', path, N_fine); 
        load (data_fine, 'valve') 

        valve_fine       = valve; 
        leaflet_fine     = valve_fine.leaflets(1); 
        X_fine           = leaflet_fine.X; 

        % final iteration 
        if (i==iterations)
            N_fine
            newton_error_difference_equations(i+1) = total_global_err(leaflet_fine); 
            
        end 
        
        % fine mesh interpolated to align with coarse mesh 
        X_fine_interp  = zeros(size(X_coarse)); 

        % for all internal points on coarse mesh 
        % find coordinates in u,v plane 
        % interpolate fine mesh to this location 
        % with bilinear interpolation 
        % only do so if all four points on the fine mesh are valid 
        %for j=1:j_max_anterior_coarse 
        for j=1:j_max_coarse
        %for j=(j_max_anterior_coarse+1):j_max_coarse 
            for k=1:k_max_coarse
                
                j_fine = 2*j; 
                k_fine = 2*k - 1; 
                
                X_fine_interp(:,j,k) = X_fine(:,j_fine,k_fine); 
                
%                 % check at the ring to make sure checker gets desired order 
%                 % piecewise linear interpolant (with exact application of boundary conditions)
%                 % should give secon order on ring 
%                 if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && (k == k_max_coarse) 
%                     
%                     % fprintf('j = %d, j_fine_below = %d, u = %f, s = %f\n', j, j_fine_below, u, s); 
%                     X_ring_coarse(:,j) = X_coarse(:,j,k); 
% 
%                     if s == 0
%                         % edge case, s==0 may be the last point (and there is no point above to be read)
%                         X_ring_interp(:,j) = X_fine(:,j_fine_below,k_max_fine); 
%                     else 
%                         X_ring_interp(:,j) = (1-s) * X_fine(:,j_fine_below,k_max_fine) + s*X_fine(:,j_fine_above,k_max_fine); 
%                     end 
%                 end   
%                 
                
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

                diff_one_components(j,k) = norm(X_coarse(:,j,k) - X_fine_interp(:,j,k),1); 
                diff_two_components(j,k) = norm(X_coarse(:,j,k) - X_fine_interp(:,j,k),2); 
                diff_inf_components(j,k) = norm(X_coarse(:,j,k) - X_fine_interp(:,j,k),'inf'); 
            end 
        end 
        
        % For results to be comparable across vectors of different sizes
        % Need to scale norms so that they are double integrals 

        differences = X_coarse(:) - X_fine_interp(:);  

        % du^2 measure for double integral 
        diff_1(i)    = du_coarse*du_coarse * norm(differences, 1); 

        % on 2 norm gets a root on the measure 
        diff_2(i)    = du_coarse           * norm(differences, 2); 

        % no weight on sup norm version 
        diff_inf(i)  =                       norm(differences, 'inf'); 
        
        
        pointwise_plots = true; 
        if pointwise_plots && (i == iterations)
            n_colors = 500; 
            extended = true; 
            colormap(make_colormap(n_colors, extended)); 

            cmap = colormap;
            n_colors = size(cmap,1); 

            max_diff_1_norm   = max(max(diff_one_components)); 
            max_diff_2_norm   = max(max(diff_two_components)); 
            max_diff_inf_norm = max(max(diff_inf_components)); 
    
            colors_1_norm     = zeros(j_max_coarse,k_max_coarse,3); 
            colors_2_norm     = zeros(j_max_coarse,k_max_coarse,3); 
            colors_inf_norm   = zeros(j_max_coarse,k_max_coarse,3); 
                        
            for j=1:j_max_coarse
                for k=1:k_max_coarse
                    
                    color_idx = floor(n_colors * diff_one_components(j,k)/max_diff_1_norm);        
                    if color_idx == 0
                        color_idx = 1; 
                    end 
                    colors_1_norm(j,k,:) = cmap(color_idx,:); 
                    
                    color_idx = floor(n_colors * diff_two_components(j,k)/max_diff_2_norm); 
                    if color_idx == 0
                        color_idx = 1; 
                    end 
                    colors_2_norm(j,k,:) = cmap(color_idx,:); 
                    
                    color_idx = floor(n_colors * diff_inf_components(j,k)/max_diff_inf_norm); 
                    if color_idx == 0
                        color_idx = 1; 
                    end
                    colors_inf_norm(j,k,:) = cmap(color_idx,:);
                    
                end 
            end 
            

            x_component = squeeze(X_coarse(1,:,:)); 
            y_component = squeeze(X_coarse(2,:,:)); 
            z_component = squeeze(X_coarse(3,:,:)); 

            width = 1.0; 
            
            figure; 
            surf(x_component, y_component, z_component, colors_1_norm, 'LineWidth',width);
            title('1 norm componentwise diffs')            
            axis equal 
            axis auto 
            hold on 
            
            colormap(make_colormap(n_colors, extended)); 
            n_ticks = 5; 
            tick_array = linspace(0,1,n_ticks); 
            tick_labels = {}; 
            for n=1:length(tick_array)
                tick=tick_array(n); 
                tension = tick * max_diff_1_norm; 
                tick_labels{n} = sprintf('%.1e', tension); 
            end 
            cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels);
            
            figure; 
            surf(x_component, y_component, z_component, colors_2_norm, 'LineWidth',width);
            title('2 norm componentwise diffs')
            colorbar
            axis equal 
            axis auto 
            hold on 
            
            colormap(make_colormap(n_colors, extended)); 
            n_ticks = 5; 
            tick_array = linspace(0,1,n_ticks); 
            tick_labels = {}; 
            for n=1:length(tick_array)
                tick=tick_array(n); 
                tension = tick * max_diff_1_norm; 
                tick_labels{n} = sprintf('%.1e', tension); 
            end 
            cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels);
            
            figure; 
            surf(x_component, y_component, z_component, colors_inf_norm, 'LineWidth',width);
            title('inf norm componentwise diffs')
            colorbar
            axis equal 
            axis auto 
            hold on 
            
            colormap(make_colormap(n_colors, extended)); 
            n_ticks = 5; 
            tick_array = linspace(0,1,n_ticks); 
            tick_labels = {}; 
            for n=1:length(tick_array)
                tick=tick_array(n); 
                tension = tick * max_diff_1_norm; 
                tick_labels{n} = sprintf('%.1e', tension); 
            end 
            cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels);

        end 
    

    end
    
%     if ring_convergence_sanity_check
%         diff_1_ring 
%         diff_2_ring
%         diff_inf_ring 
% 
% 
%         fprintf('Comparisons of N,2N and 2N,4N points\n'); 
%         N = N_values(1); 
%         for i = 1:iterations-1
% 
%             order_1    = diff_1_ring(i) / diff_1_ring(i+1); 
%             order_2    = diff_2_ring(i) / diff_2_ring(i+1);  
%             order_inf  = diff_inf_ring(i) / diff_inf_ring(i+1);  
% 
%             fprintf('%d,%d,%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N, 2*N, 4*N, order_1, order_2, order_inf);
% 
%             N = 2*N;
%         end
%     end 
    

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
    
end 