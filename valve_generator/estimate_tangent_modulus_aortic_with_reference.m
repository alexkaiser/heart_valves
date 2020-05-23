function [sigma_circ, sigma_rad, sigma_circ_mean, sigma_rad_mean]  = estimate_tangent_modulus_aortic_with_reference(leaflet, thickness, plots)
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
    
    sigma_circ = zeros(j_max, k_max); 
    sigma_rad  = zeros(j_max, k_max); 
        
    collagen_constitutive_circ = leaflet.collagen_constitutive_circ; 
    collagen_constitutive_rad  = leaflet.collagen_constitutive_rad; 
    
    if ~exist('plots', 'var')
        plots = false; 
    end 
    
    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k)
                
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

                        n_valid_cross_circ = 0; 
                        
                        % cross wise direction springs 
                        for k_nbr_tmp_cross = [k-1,k+1]
                            [valid_cross j_nbr_cross k_nbr_cross] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp_cross);
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
                        sigma_circ(j,k) = sigma_circ(j,k) + tension_grad / (len_element_u_type_temp * thickness);                         

                    elseif valid && target_spring 
                        error('No j direction targets allowed'); 
                    end 

                end  
                
                % average forces collected 
                if element_count_circ > 0
                    sigma_circ(j,k) = sigma_circ(j,k)/element_count_circ; 
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

                        n_valid_cross_rad = 0; 
                        
                        % cross wise direction springs 
                        for j_nbr_tmp_cross = [j-1,j+1]
                            [valid_cross j_nbr_cross k_nbr_cross] = get_indices(leaflet, j, k, j_nbr_tmp_cross, k_nbr_tmp);
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
                        sigma_rad(j,k) = sigma_rad(j,k) + tension_grad / (len_element_v_type_temp * thickness);
                    
                    end

                end                

                % average forces collected 
                if element_count_rad > 0
                    sigma_rad(j,k) = sigma_rad(j,k)/element_count_rad; 
                end 

            end
        end  
    end
    
    sigma_circ_linear = sigma_circ(:); 
    sigma_circ_valid  = sigma_circ_linear(find(sigma_circ_linear ~= 0)); 
    
    sigma_rad_linear = sigma_rad(:); 
    sigma_rad_valid  = sigma_rad_linear(find(sigma_rad_linear ~= 0)); 
    
    sigma_circ_mean = mean(sigma_circ_valid); 
    sigma_rad_mean  = mean(sigma_rad_valid); 

    if plots
        for circ_temp = [true, false]
            
            if circ_temp
                sigma = sigma_circ; 
            else
                sigma = sigma_rad; 
            end 
            figure; 
            
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

            max_sigma = max(max(sigma));

            tick_max = max_sigma; 
            colors = zeros(j_max,k_max,3); 
            for j=1:j_max
                for k=1:k_max
                    color_idx = floor(n_colors * sigma(j,k) / max_sigma); 
                    if color_idx == 0
                        color_idx = 1; 
                    end 
                    if color_idx > n_colors
                        color_idx = n_colors; 
                    end         
                    colors(j,k,:) = cmap(color_idx,:); 
                end 
            end 

            x_component = squeeze(X_copy(1,:,:)); 
            y_component = squeeze(X_copy(2,:,:)); 
            z_component = squeeze(X_copy(3,:,:)); 
            % surf(x_component, y_component, z_component, colors, 'edgecolor', 'none');
            surf(x_component, y_component, z_component, colors);
            % colormap(make_colormap(n_colors, extended)); 
            
            n_ticks = 5; 
            tick_array = linspace(0,1,n_ticks); 
            tick_labels = {}; 
            for n=1:length(tick_array)
                tick=tick_array(n); 
                tmp = tick * tick_max; 
                tick_labels{n} = sprintf('%.1e', tmp); 
            end 
            colorbar('Ticks', tick_array, 'TickLabels', tick_labels);

            view(14,14)
            axis equal 
            
            if circ_temp
                title('Tangent modulus at current strain, circumferential'); 
            else 
                title('Tangent modulus at current strain, radial'); 
            end 
            
        end 
    end
    
    
end

 
