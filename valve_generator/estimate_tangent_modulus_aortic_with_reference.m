function [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean, fig, stress_circ_mean, stress_rad_mean] ...
    = estimate_tangent_modulus_aortic_with_reference(leaflet, thickness, fig, fiber_stride, stride_offset_j, circ, rad, ratio, max_plot_cap, plot_stress)
    % estimates the tangent modulus for current strain 
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

    X_current          = leaflet.X; 
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    is_internal        = leaflet.is_internal; 
    is_bc              = leaflet.is_bc; 
    R_u                = leaflet.R_u;
    k_u                = leaflet.k_u;
    R_v                = leaflet.R_v;
    k_v                = leaflet.k_v;
    N_each             = leaflet.N_each; 
    
    if isfield(leaflet, 'N_leaflets')
        N_leaflets         = leaflet.N_leaflets; 
    else 
        % trileaflet default
        N_leaflets         = 3; 
    end 
    % tangent modulus, slope of stress 
    sigma_circ = zeros(j_max, k_max); 
    sigma_rad  = zeros(j_max, k_max); 
        
    % values of stress 
    stress_circ = zeros(j_max, k_max); 
    stress_rad  = zeros(j_max, k_max); 
    
    collagen_constitutive_circ = leaflet.collagen_constitutive_circ; 
    collagen_constitutive_rad  = leaflet.collagen_constitutive_rad; 
    
    if ~exist('fig', 'var')
        fig = figure;  
    end 

    if ~exist('fiber_stride', 'var')
        fprintf('Using default fiber stride of 1')
        fiber_stride = 1; 
    end 

    if ~exist('stride_offset_j', 'var')
        stride_offset_j = 0; 
    end 

    hold on;

    if ~exist('circ', 'var')
        circ = false; 
    end 
    
    if ~exist('rad', 'var')
        rad = false; 
    end 
    
    if ~exist('ratio', 'var')
        ratio = false; 
    end 
    
    plots = (circ || rad || ratio); 
    
    if (circ && rad) || (circ && ratio) || (rad && ratio)
        error('too many plots requested at once')
    end 
    
    if ~exist('max_plot_cap', 'var')
        max_plot_cap = inf; 
    end 
    
    if rad && isinf(max_plot_cap)
        radial_autoscale = true; 
    end 
    
    if ~exist('plot_stress', 'var')
        plot_stress = false;
    end 
    
    one_leaflet = true; 
    if one_leaflet 
        j_min_plot = N_each; 
        j_max_plot = 2*N_each; 
    end 
    
    
    
    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            
            if (mod(j,fiber_stride) == (1 + stride_offset_j)) || (fiber_stride == 1)
                output_tmp_k = true; 
            else 
                output_tmp_k = false; 
            end 

            if (mod(k,fiber_stride) == 1) || (fiber_stride == 1)
                output_tmp_j = true; 
            else 
                output_tmp_j = false; 
            end     

            if one_leaflet 
                if (j < j_min_plot) || (j > j_max_plot)
                    output_tmp_j = false; 
                    output_tmp_k = false; 
                end 
            end 
                
            X = X_current(:,j,k); 

            element_count_circ = 0; 

            % u type fibers 
            for j_nbr_tmp = [j-1,j+1]

                len_element_u_type_temp = 0; 

                k_nbr_tmp = k; 

                [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid && (~target_spring) && (~target_k_no_j_spring)
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    tension_grad = tension_derivative_with_reference_wrt_strain(X, X_nbr, R_u(j_spr,k_spr), k_u(j_spr,k_spr), leaflet, collagen_constitutive_circ); 
                    tension = tension_with_reference(X, X_nbr, R_u(j_spr,k_spr), k_u(j_spr,k_spr), leaflet, collagen_constitutive_circ); 

                    n_valid_cross_circ = 0; 

                    % cross wise direction springs 
                    for k_nbr_tmp_cross = [k-1,k+1]
                        [valid_cross j_nbr_cross k_nbr_cross] = get_indices(leaflet, j, k, j, k_nbr_tmp_cross);
                        if valid_cross
                            n_valid_cross_circ = n_valid_cross_circ + 1; 
                            X_nbr_cross = X_current(:,j_nbr_cross,k_nbr_cross); 
                            len_element_u_type_temp = len_element_u_type_temp + 0.5 * norm(X - X_nbr_cross); 
                        end 
                    end 

                    if n_valid_cross_circ == 0 
                        error('found spring without any valid neighbors in cross direction')
                    end 

                    element_count_circ = element_count_circ + 1; 
                    if ~isnan(tension_grad)
                        sigma_circ(j,k) = sigma_circ(j,k) + tension_grad / (len_element_u_type_temp * thickness);
                    end 
                    
                    if ~isnan(tension)
                        stress_circ(j,k) = stress_circ(j,k) + tension / (len_element_u_type_temp * thickness);
                    end 
                    
                    if one_leaflet 
                        if (j_nbr < j_min_plot) || (j_nbr > j_max_plot)
                            output_tmp_j = false; 
                            output_tmp_k = false; 
                        end 
                    end 
                    
                    x_vals = [X(1), X_nbr(1)]; 
                    y_vals = [X(2), X_nbr(2)]; 
                    z_vals = [X(3), X_nbr(3)]; 
                    if output_tmp_j && (circ || ratio) && plots 
                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 

                elseif valid && target_spring 
                    error('No j direction targets allowed'); 
                end 

            end  

            % average forces collected 
            if element_count_circ > 0
                sigma_circ(j,k)  = sigma_circ(j,k)/element_count_circ; 
                stress_circ(j,k) = stress_circ(j,k)/element_count_circ; 
            end 

            element_count_rad = 0; 

            % v type fibers 
            for k_nbr_tmp = [k-1,k+1]

                len_element_v_type_temp = 0; 

                j_nbr_tmp = j; 

                [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid && (~target_spring)
                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    tension_grad = tension_derivative_with_reference_wrt_strain(X, X_nbr, R_v(j_spr,k_spr), k_v(j_spr,k_spr), leaflet, collagen_constitutive_rad); 
                    tension = tension_with_reference(X, X_nbr, R_v(j_spr,k_spr), k_v(j_spr,k_spr), leaflet, collagen_constitutive_rad); 

                    n_valid_cross_rad = 0; 

                    % cross wise direction springs 
                    for j_nbr_tmp_cross = [j-1,j+1]
                        [valid_cross j_nbr_cross k_nbr_cross] = get_indices(leaflet, j, k, j_nbr_tmp_cross, k);
                        if valid_cross
                            n_valid_cross_rad = n_valid_cross_rad + 1; 
                            X_nbr_cross = X_current(:,j_nbr_cross,k_nbr_cross); 
                            len_element_v_type_temp = len_element_v_type_temp + 0.5 * norm(X - X_nbr_cross); 
                        end 
                    end 

                    if n_valid_cross_rad == 0 
                        error('found spring without any valid neighbors in cross direction')
                    end 

                    element_count_rad = element_count_rad + 1; 
                    if ~isnan(tension_grad)
                        sigma_rad(j,k) = sigma_rad(j,k) + tension_grad / (len_element_v_type_temp * thickness);
                    end 
                    if ~isnan(tension)
                        stress_rad(j,k) = stress_rad(j,k) + tension / (len_element_v_type_temp * thickness);
                    end 
                    
                    x_vals = [X(1), X_nbr(1)]; 
                    y_vals = [X(2), X_nbr(2)]; 
                    z_vals = [X(3), X_nbr(3)]; 
                    % plot local fiber if included 
                    if output_tmp_k && (rad || ratio) && plots  
                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 

                end

            end                

            % average forces collected 
            if element_count_rad > 0
                sigma_rad(j,k)  = sigma_rad(j,k)/element_count_rad; 
                stress_rad(j,k) = stress_rad(j,k)/element_count_rad; 
            end 

            
        end  
    end
    
    sigma_circ_linear = sigma_circ(:); 
    sigma_circ_valid  = sigma_circ_linear(find(sigma_circ_linear ~= 0)); 
    
    sigma_rad_linear = sigma_rad(:); 
    sigma_rad_valid  = sigma_rad_linear(find(sigma_rad_linear ~= 0)); 
    
    sigma_circ_mean = mean(sigma_circ_valid); 
    sigma_rad_mean  = mean(sigma_rad_valid); 

    stress_circ_linear = stress_circ(:); 
    stress_circ_valid  = stress_circ_linear(find(stress_circ_linear ~= 0)); 
    
    stress_rad_linear = stress_rad(:); 
    stress_rad_valid  = stress_rad_linear(find(stress_rad_linear ~= 0)); 
    
    stress_circ_mean = mean(stress_circ_valid); 
    stress_rad_mean  = mean(stress_rad_valid); 
    
    if plots
        
        if plot_stress
            if circ
                sigma = stress_circ; 
            elseif rad 
                sigma = stress_rad; 
            elseif ratio 
                sigma = stress_circ ./ stress_rad; 
            end 
        else 
            % tangent mod default 
            if circ
                sigma = sigma_circ; 
            elseif rad 
                sigma = sigma_rad; 
            elseif ratio 
                sigma = sigma_circ ./ sigma_rad; 
            end 
        end 
        
        if ratio 
            mean_ratio = 0; 
            valid_ratio = 0; 
            
            for j=1:j_max
                for k=1:k_max
                    if ~((sigma(j,k) == 0) || isnan(sigma(j,k)) || isinf(sigma(j,k)))
                        mean_ratio = mean_ratio + sigma(j,k); 
                        valid_ratio = valid_ratio + 1; 
                    end 
                end 
            end 
            
            mean_ratio = mean_ratio / valid_ratio;             
            fprintf('mean ratio, ratio first = %f', mean_ratio); 
        end 
        
        n_colors = 500;
        extended = true; 
        colormap(make_colormap(n_colors, extended)); 
        cmap = colormap;
        n_colors = size(cmap,1); 

        % plot the actual surface 
        X_copy      = leaflet.X; 

        % NaN mask in the copy 
        for j=1:j_max
            for k=1:k_max
                if ~(is_internal(j,k) || is_bc(j,k))
                   X_copy(:,j,k) = NaN;  
                end
            end 
        end

        if any(any(any(isnan(X_copy))))
            error('all points should be internal or bcs...')
        end 

        if one_leaflet 

            X_copy(:,1:(j_min_plot-1),:) = NaN;    
            X_copy(:,(j_max_plot+1):end,:) = NaN;

            outline_cleanup = true; 
            if outline_cleanup 

                % horizontal at bottom and free edge 
                for k=[1,k_max]
                    for j=j_min_plot:(j_max_plot-1)

                        X = X_copy(:,j,k); 
                        X_nbr = X_copy(:,j+1,k); 

                        x_vals = [X(1), X_nbr(1)]; 
                        y_vals = [X(2), X_nbr(2)]; 
                        z_vals = [X(3), X_nbr(3)]; 

                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 
                end 

                % vertical at commissures 
                for k=1:(k_max-1)
                    for j=[j_min_plot,j_max_plot]

                        X = X_copy(:,j,k); 
                        X_nbr = X_copy(:,j,k+1); 

                        x_vals = [X(1), X_nbr(1)]; 
                        y_vals = [X(2), X_nbr(2)]; 
                        z_vals = [X(3), X_nbr(3)]; 

                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 
                end 
            end 
        end 
        
        if isinf(max_plot_cap)
            max_plot_cap = max(max(sigma));
        end 
            
        % at commissures in vertical direction
        % where radial tension is not defined 
        % take the mean of neighbors 
        comm_color_patch = true; 
        if (rad || ratio) && comm_color_patch
            
            for j=(N_each * (1:N_leaflets))
                for k=1:k_max
                    
                    if ~((sigma(j,k) == 0) || isnan(sigma(j,k)) || isinf(sigma(j,k)))
                        error('trying to patch valid color'); 
                    end 
                                            
                    for j_nbr_tmp = [j-1,j+1]
                        k_nbr_tmp = k; 
                        [valid j_nbr k_nbr j_spr k_spr target_spring target_k_no_j_spring] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                        
                        if ~valid 
                            error('should always have two neighbors here'); 
                        end 
                        
                        sigma(j,k) = sigma(j,k) + sigma(j_nbr,k_nbr); 
                        
                    end 
                    
                    sigma(j,k) = sigma(j,k) / 2; 
                    
                end 
            end 

        end 
        
        tick_max = max_plot_cap; 
        colors = zeros(j_max,k_max,3); 
        for j=1:j_max
            for k=1:k_max
                color_idx = floor(n_colors * sigma(j,k) / max_plot_cap); 
                if (color_idx == 0) || (isnan(color_idx))
                    color_idx = 1; 
                end 
                if color_idx > n_colors
                    color_idx = n_colors; 
                end         
                colors(j,k,:) = cmap(color_idx,:); 
            end 
        end 
        
        circ_off_by_one_adjust = true; 
        if (circ || ratio) && circ_off_by_one_adjust 
            % colors are read from the minimum j,k index for the square 
            % so shift down by one 
            colors_temp = zeros(size(colors)); 

            for j=1:j_max 
                for k=1:(k_max-1)
                    colors_temp(j,k,:) = colors(j,k+1,:); 
                end 
            end 

            colors = colors_temp; 
        end 

        x_component = squeeze(X_copy(1,:,:)); 
        y_component = squeeze(X_copy(2,:,:)); 
        z_component = squeeze(X_copy(3,:,:)); 
        surf(x_component, y_component, z_component, colors, 'edgecolor', 'none');
        % surf(x_component, y_component, z_component, colors);
        % colormap(make_colormap(n_colors, extended)); 

        n_ticks = 5; 
        tick_array = linspace(0,1,n_ticks); 
        tick_labels = {}; 
        for n=1:length(tick_array)
            tick=tick_array(n); 
            tmp = tick * tick_max; 
            if ratio 
                tick_labels{n} = sprintf('%.0f', tmp); 
            else 
%                 exponent_removed = 10^floor(log10(tick_max)); 
%                 tick_labels{n} = sprintf('%.1f', tmp/exponent_removed); 
                tick_labels{n} = sprintf('%.1e', tmp); 
            end 
        end 
        
        colorbar_on = true; 
        if colorbar_on
            colorbar('Ticks', tick_array, 'TickLabels', tick_labels);
        end 
        
        colorbar_figure = false; 
        if colorbar_figure 
            fig_colorbar = figure; 

            colormap(make_colormap(n_colors, extended)); 

            cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels); 

            fontsize = 24; 
            ax = gca; 
            ax.FontSize = fontsize;
            cbar.Label.FontSize = fontsize; 
            cbar.Label.Rotation = 0;
            cbar.Label.Position = [0.4 1.2];

            grid off 
            axis off 

            if circ 
                bar_name = 'colorbar_only_circ_tangent_mod'; 
            elseif rad && exist('radial_autoscale', 'var') && radial_autoscale
                exponent_removed 
                bar_name = 'colorbar_only_radial_tangent_mod_autoscale'; 
            elseif rad 
                bar_name = 'colorbar_only_radial_tangent_mod'; 
            elseif ratio 
                bar_name = 'colorbar_only_ratio_tangent_mod'; 
            else
                error('incompatible format arguments')
            end 

            print(fig_colorbar, '-depsc', bar_name); 
            close(fig_colorbar);    
            % reset current figure 
            figure(fig); 
        end 

        axis equal 
        axis off 
        axis tight

    end 
end


 
