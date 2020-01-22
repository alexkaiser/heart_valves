function [fig] = total_tension_surf_plot_aortic(leaflet, fiber_output, fiber_stride, stride_offset_j, circ, rad, ratio, height_plot, fig)
% 
% Plots leaflet with fibers 
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


% 
% Note following example: 
% 
% https://www.mathworks.com/matlabcentral/answers/34750-plot3-color
% 
% figure
% cmap = colormap;
% % change c into an index into the colormap
% % min(c) -> 1, max(c) -> number of colors
% c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
% % make a blank plot
% plot3(x,y,z,'linestyle','none')
% % add line segments
% for k = 1:(length(x)-1)
%     line(x(k:k+1),y(k:k+1),z(k:k+1),'color',cmap(c(k),:))
% end
% colorbar
% 
% caxis([ min(c) , max(c)]) % colorbar limits 



% basic idea -- copy difference euqation code
% every time there is a link and the tension is computed, plot that link 
% with the color as tension 
% 

% set the max tension on the color bar to be the max of the alpha or beta 
% as appropriate for this family of fibers 

% consider also the chordae... 
% this have no notion of force density and must be scaled to be compatible
% with the tensions in the leaflet 



% 
% Evaluation of the global difference equations at j,k
% Requires reference configuration R 
% 
% Input
%     leaflet    Current parameters 
%
% Output
%     F        Values of all difference equation, 3 by triangular array 
% 

X_current              = leaflet.X; 
alpha                  = leaflet.alpha; 
beta                   = leaflet.beta; 
c_dec_radial           = leaflet.c_dec_radial; 
c_dec_circumferential  = leaflet.c_dec_circumferential; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 

if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
    decreasing_tension    = true;  
else 
    decreasing_tension    = false;  
end 

if isfield(leaflet, 'tension_debug') && leaflet.tension_debug
    tension_debug = true; 
else 
    tension_debug = false; 
end 

if ~exist('fig', 'var')
    fig = figure;  
end 

if ~exist('ratio', 'var')
    ratio = false; 
end 

if ~exist('height_plot', 'var')
    height_plot = false; 
end 


% quick argument checking 
if circ && rad
elseif circ  
elseif rad 
elseif ratio
else 
    error('invalid arguments');    
end

if exist('fiber_output', 'var') && fiber_output    
    if ~exist('fiber_stride', 'var')
        fprintf('Using default fiber stride of 1')
        fiber_stride = 1; 
    end 
    
    if ~exist('stride_offset_j', 'var')
        stride_offset_j = 0; 
    end 
    
end 

hold on; 

n_colors = 500; 
extended = true; 
colormap(make_colormap(n_colors, extended)); 


cmap = colormap;
n_colors = size(cmap,1); 

max_tension_circ = du * max(alpha(:)); 
max_tension_radial = du * max(beta(:)); 
if circ && rad 
    max_tension = max_tension_circ + max_tension_radial; % 1.7 * 0.7 * max(max_tension_circ, max_tension_radial); 
else 
    max_tension = max(max_tension_circ, max_tension_radial); 
end
tick_max = max_tension; 


% color leaflet by total of tension 
% required to be m,n,3 (in annoying contrast to my normal convention)
colors   = nan * zeros(j_max,k_max,3); 
% counts the number of tensions added to any given node 
tension_circ  =  zeros(j_max,k_max); 
num_nbrs_circ =  zeros(j_max,k_max); 
tension_rad   =  zeros(j_max,k_max); 
num_nbrs_rad  =  zeros(j_max,k_max); 

ratio_ptwise = zeros(j_max,k_max); 

% Internal leaflet part 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) 
            
            if fiber_output && ((mod(j,fiber_stride) == 1) || (fiber_stride == 1))
                output_tmp_k = true; 
            else 
                output_tmp_k = false; 
            end 

            if fiber_output && ((mod(k,fiber_stride) == (1 + stride_offset_j)) || (fiber_stride == 1))
                output_tmp_j = true; 
            else 
                output_tmp_j = false; 
            end     
            

            X = X_current(:,j,k); 

            F_tmp = zeros(3,1); 


            % u type fibers 
            for j_nbr_tmp = [j-1,j+1]

                k_nbr_tmp = k; 

                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid 

                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    alpha_tmp     = alpha(j_spr,k_spr); 
                    c_dec_tension = c_dec_circumferential(j_spr,k_spr); 

                    tension = alpha_tmp;  

                    if decreasing_tension && (alpha_tmp ~= 0)
                        tension = tension + alpha_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                    end 

                    if tension_debug
                        dec = tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                        fprintf('tension = %e, dec_tension = %f, (j,k) = (%d, %d) circ\n', tension, dec, j, k); 
                    end 

                    % ensure that the circumferential figure is current 
                    % figure(fig); 
                    
                    x_vals = [X(1), X_nbr(1)]; 
                    y_vals = [X(2), X_nbr(2)]; 
                    z_vals = [X(3), X_nbr(3)]; 
                    
                    % fraction of maximum tension gives fraction of way
                    % through color bar
                    tension_circ(j,k)  = tension_circ(j,k) + du * tension; 
                    num_nbrs_circ(j,k) = num_nbrs_circ(j,k) + 1; 
                    
                    if output_tmp_j && circ 
                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 

                    F_tmp = F_tmp + du * tension * (X_nbr-X)/norm(X_nbr-X); 

                end 
            end 

            % v type fibers 
            for k_nbr_tmp = [k-1,k+1]

                j_nbr_tmp = j; 

                [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid

                    X_nbr = X_current(:,j_nbr,k_nbr); 

                    beta_tmp      = beta(j_spr,k_spr); 
                    c_dec_tension = c_dec_radial(j_spr,k_spr); 

                    tension = beta_tmp; 

                    if decreasing_tension && (beta_tmp ~= 0)
                        tension = tension + beta_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                    end

                    if tension_debug
                        dec = tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                        fprintf('tension = %e, dec_tension = %f, (j,k) = (%d, %d) radial\n', tension, dec, j, k); 
                    end 

                    % ensure that the circumferential figure is current 
                    % figure(fig); 
                    
                    x_vals = [X(1), X_nbr(1)]; 
                    y_vals = [X(2), X_nbr(2)]; 
                    z_vals = [X(3), X_nbr(3)]; 
                    
                    % fraction of maximum tension gives fraction of way
                    % through color bar 
                    tension_rad(j,k) = tension_rad(j,k) + du * tension; 
                    num_nbrs_rad(j,k) = num_nbrs_rad(j,k) + 1;  

                    % plot local fiber if included 
                    if output_tmp_k && rad 
                        plot3(x_vals,y_vals,z_vals,'k'); 
                    end 
                    
                end 

            end 

            
            % update color map for sum of local tensions 
            % tensions either take the exact value (if they only are connected to one node)
            % or the average on the connecting segments 
            tension_sum = 0; 
            
            if num_nbrs_circ(j,k) == 1
                circ_temp =       tension_circ(j,k); 
            elseif num_nbrs_circ(j,k) == 2
                circ_temp = 0.5 * tension_circ(j,k); 
            else 
                error('Improper number of circumferential tension nodes'); 
            end    

            if num_nbrs_rad(j,k) == 1
                rad_temp =       tension_rad(j,k); 
            elseif num_nbrs_rad(j,k) == 2
                rad_temp = 0.5 * tension_rad(j,k); 
            else 
                error('Improper number of radial tension nodes'); 
            end    
            
            ratio_ptwise(j,k) = circ_temp / rad_temp; 
            
            if circ 
                tension_sum = tension_sum + circ_temp; 
            end 
            if rad 
                tension_sum = tension_sum + rad_temp; 
            end 
            
            tension_frac_of_max = tension_sum/max_tension; 
            color_idx = ceil(tension_frac_of_max * n_colors);

            if color_idx == 0
                color_idx = 1; 
            end 

            if color_idx > n_colors
                %warning('max color at free edge')
                color_idx = n_colors; 
            end 
            
            colors(j,k,:) = cmap(color_idx,:); 

        end
    end 
end

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

if ratio 
    max_ratio = max(max(ratio_ptwise));
    tick_max = max_ratio; 
    colors = zeros(j_max,k_max,3); 
    for j=1:j_max
        for k=1:k_max
            color_idx = floor(n_colors * ratio_ptwise(j,k) / max_ratio); 
            if color_idx == 0
                color_idx = 1; 
            end 
            if color_idx > n_colors
                color_idx = n_colors; 
            end         
            colors(j,k,:) = cmap(color_idx,:); 
        end 
    end 
end 


if height_plot
    if ~ratio
        error('height only for ratio')
    end
    
    colormap default; 
    
    height_cap = 10; 
    ratio_ptwise = min(height_cap, ratio_ptwise); 
    
    % patch bottom, copy from one up 
    ratio_ptwise(:,1) = nan; 
    
    x_mesh = du * (1:j_max); 
    y_mesh = du * (1:k_max); 
    [x_component y_component] = meshgrid(x_mesh, y_mesh); 
    x_component = x_component'; 
    y_component = y_component';     
    colorbar;
    surf(x_component, y_component, ratio_ptwise)
    axis([0 max(x_mesh) 0 max(y_mesh) 0 10])
    view(30,30); 
     
else
    
    n_ticks = 5; 
    tick_array = linspace(0,1,n_ticks); 
    tick_labels = {}; 
    for n=1:length(tick_array)
        tick=tick_array(n); 
        tmp = tick * tick_max; 
        tick_labels{n} = sprintf('%.1e', tmp); 
    end 
    colorbar('Ticks', tick_array, 'TickLabels', tick_labels);
    
    x_component = squeeze(X_copy(1,:,:)); 
    y_component = squeeze(X_copy(2,:,:)); 
    z_component = squeeze(X_copy(3,:,:)); 
    % surf(x_component, y_component, z_component, colors, 'edgecolor', 'none');
    surf(x_component, y_component, z_component, colors);
    % colormap(make_colormap(n_colors, extended)); 
end 
    
    
if circ && rad
    title('total tension')
elseif circ  
    title('circ tension')
elseif rad 
    title('radial tension')
elseif ratio
    title(sprintf('ratio circ/rad tension, max ratio = %f', max_ratio)); 
else 
    error('invalid arguments')    
end


