function [fig] = total_tension_surf_plot(leaflet, anterior, fiber_output, fiber_stride, stride_offset_j, circ, rad, fig)
% 
% Plots leaflet with fibers 
% 
% 


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
p_0                    = leaflet.p_0; 
alpha                  = leaflet.alpha; 
beta                   = leaflet.beta; 
c_dec_radial           = leaflet.c_dec_radial; 
c_dec_circumferential  = leaflet.c_dec_circumferential; 
chordae                = leaflet.chordae; 
chordae_idx            = leaflet.chordae_idx; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
num_trees              = leaflet.num_trees; 

% allows for anterior/posterior only 
j_range_anterior   = leaflet.j_range_anterior; 
j_range_right_comm = leaflet.j_range_right_comm; 
j_range_posterior  = leaflet.j_range_posterior; 
j_range_left_comm  = leaflet.j_range_left_comm; 

% two anterior trees for now, gross hack 
n_trees_anterior = 2; 
if anterior 
    j_range    = j_range_anterior; 
    tree_range = 1:n_trees_anterior; 
else 
    j_range = [j_range_right_comm, j_range_posterior, j_range_left_comm]; 
    tree_range = (n_trees_anterior+1):num_trees; 
end 


if isfield(leaflet, 'periodic_j')
    periodic_j = leaflet.periodic_j; 
else
    periodic_j = zeros(k_max,1); 
end 

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

if ~exist('circ', 'var')
    circ = true;  
end 

if ~exist('rad', 'var')
    rad = true;  
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

if (~circ) && (~rad)
    error('Have to plot at least one fiber family for this function to do anything')
end 

% x_tmp = X_current(1,:,:)
plot3(0,0,0,'linestyle','none'); 
hold on; 

% inverted hot colormap 
% c = hot;
% c = flipud(c);
% colormap(c);

% wide color range, red at top 
% colormap(jet); 

% range = 100; 
% jet_wide_range = jet(range); 
% min_idx = 30; 
% jet_cropped = jet_wide_range(min_idx:end,:); 
% colormap(jet_cropped)

% for i=1:(range/4)
%      jet_wide_range(i,:) = (1 - (i-1)/range) * jet_wide_range(i,:)/norm(jet_wide_range(i,:)); 
% end 
% colormap(jet_wide_range); 

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

% doesn't make any difference here
% max_tension_sum = du * max(alpha(:) + beta(:)); 

% cbar = colorbar; 

n_ticks = 5; 
tick_array = linspace(0,1,n_ticks); 
tick_labels = {}; 
for i=1:length(tick_array)
    tick=tick_array(i); 
    tension = tick * max_tension * 1e-3; 
    tick_labels{i} = sprintf('%.1f', tension); 
end 

% cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels); 
% cbar.Label.String = 'Tension (K Dyne)'; 

colorbar_figure = true; 
if colorbar_figure 
    fig_colorbar = figure; 
    
    colormap(make_colormap(n_colors, extended)); 
    

    cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels); 
    
    cbar.Label.String = {'Tension','(\cdot 10^3 dynes)'};
    
    fontsize = 16; 
    ax = gca; 
    ax.FontSize = fontsize;
    cbar.Label.FontSize = fontsize; 
    cbar.Label.Rotation = 0;
    cbar.Label.Position = [0.4 1.2];
   
    grid off 
    axis off 
    if circ && rad 
        printfig(fig_colorbar, 'colorbar_only_total_tension'); 
    else
        printfig(fig_colorbar, 'colorbar_only_surf_one_family_tension'); 
    end
    close(fig_colorbar);    
    % reset current figure 
    figure(fig); 
end 


% color leaflet by total of tension 
% required to be m,n,3 (in annoying contrast to my normal convention)
colors   = nan * zeros(j_max,k_max,3); 
% counts the number of tensions added to any given node 
tension_circ  =  zeros(j_max,k_max); 
num_nbrs_circ =  zeros(j_max,k_max); 
tension_rad   =  zeros(j_max,k_max); 
num_nbrs_rad  =  zeros(j_max,k_max); 

% Internal leaflet part 
for j=j_range
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
%                     tension_frac_of_max = du * tension/max_tension; 
%                     color_idx = ceil(tension_frac_of_max * n_colors); 
%                     
%                     if color_idx == 0
%                         color_idx = 1; 
%                     end 
%                     
%                     if color_idx > n_colors
%                         %warning('max color, something off in color indexing')
%                         color_idx = n_colors; 
%                     end
%                     
%                     colors(j,k,:) = colors(j,k,:) + cmap(color_idx,:); 
%                     
%                     if circ
%                         plot3(x_vals,y_vals,z_vals,'color',cmap(color_idx,:)); 
%                     end 
                    
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
%                     tension_frac_of_max = du * tension/max_tension; 
%                     color_idx = ceil(tension_frac_of_max * n_colors); 
%                     
%                     if color_idx == 0
%                         color_idx = 1; 
%                     end 
%                     
%                     if color_idx > n_colors
%                         %warning('max color, something off in color indexing')
%                         color_idx = n_colors; 
%                     end 
%                     
%                     if radial 
%                         plot3(x_vals,y_vals,z_vals,'color',cmap(color_idx,:)); 
%                     end 

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

            if circ 
                if num_nbrs_circ(j,k) == 1
                    tension_sum = tension_sum + tension_circ(j,k); 
                elseif num_nbrs_circ(j,k) == 2
                    tension_sum = tension_sum + 0.5 * tension_circ(j,k); 
                else 
                    error('Improper number of circumferential tension nodes'); 
                end    
            end 
            
            if rad 
                if num_nbrs_rad(j,k) == 1
                    tension_sum = tension_sum +       tension_rad(j,k); 
                elseif num_nbrs_rad(j,k) == 2
                    tension_sum = tension_sum + 0.5 * tension_rad(j,k); 
                else 
                    error('Improper number of radial tension nodes'); 
                end    
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
            

            % current node has a chordae connection
            if chordae_idx(j,k).tree_idx

                tree_idx = chordae_idx(j,k).tree_idx; 

                [m N_chordae] = size(chordae(tree_idx).C);
                c_dec_tension_chordae = chordae(tree_idx).c_dec_chordae_leaf; 
                du_chordae = 1; 

                kappa = chordae(tree_idx).k_0;

                % index that free edge would have if on tree
                % remember that leaves are only in the leaflet
                leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                % then take the parent index of that number in chordae variables
                idx_chordae = floor(leaf_idx/2);

                X_nbr = chordae(tree_idx).C(:,idx_chordae);
                tension = kappa;  

                if decreasing_tension && (kappa ~= 0)
                    tension = tension + kappa * tension_decreasing(X, X_nbr, du_chordae, c_dec_tension_chordae); 
                end

                if tension_debug
                    dec = tension_decreasing(X, X_nbr, du, c_dec_tension_chordae); 
                    fprintf('tension = %e, dec_tension = %f, (j,k) = (%d, %d) free edge\n', tension, dec, j, k); 
                end 

                % ensure that the circumferential figure is current 
                % figure(fig); 

                x_vals = [X(1), X_nbr(1)]; 
                y_vals = [X(2), X_nbr(2)]; 
                z_vals = [X(3), X_nbr(3)]; 

                % fraction of maximum tension gives fraction of way
                % through color bar 
                tension_frac_of_max = tension/max_tension; 
                color_idx = ceil(tension_frac_of_max * n_colors);

                if color_idx == 0
                    color_idx = 1; 
                end 

                if color_idx > n_colors
                    %warning('max color at free edge')
                    color_idx = n_colors; 
                end 

                plot3(x_vals,y_vals,z_vals,'color',cmap(color_idx,:)); 
                
                F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

            end 

        end
    end 
end






% plot the actual surface 
X_copy      = leaflet.X; 
% j_max       = leaflet.j_max; 
% k_max       = leaflet.k_max; 
% is_internal = leaflet.is_internal; 
% is_bc       = leaflet.is_bc; 
% 
% j_range_anterior   = leaflet.j_range_anterior; 
% j_range_right_comm = leaflet.j_range_right_comm; 
% j_range_posterior  = leaflet.j_range_posterior; 
% j_range_left_comm  = leaflet.j_range_left_comm; 

% NaN mask in the copy 
for j=1:j_max
    for k=1:k_max
        if ~(is_internal(j,k) || is_bc(j,k))
           X_copy(:,j,k) = NaN;  
        end
    end 
end
        
x_component = squeeze(X_copy(1,j_range,:)); 
y_component = squeeze(X_copy(2,j_range,:)); 
z_component = squeeze(X_copy(3,j_range,:)); 
colors_local = colors(j_range,:,:); 
width = 1.0; 
surf(x_component, y_component, z_component, colors_local, 'edgecolor', 'none');

    

% chordae internal terms 
for tree_idx = tree_range

    C = chordae(tree_idx).C; 
    [m N_chordae] = size(C);         
    % c_dec_tension_chordae = chordae(tree_idx).c_dec_tension_chordae; 
    F_chordae(tree_idx).C = zeros(size(C));  

    % normalize this, no mesh parameters in chordae computations 
    du_chordae = 1; 

    for i=1:N_chordae

        left   = 2*i; 
        right  = 2*i + 1;
        parent = floor(i/2); 

        for nbr_idx = [left,right,parent]

            % get the neighbors coordinates, reference coordinate and spring constants
            [nbr R_nbr k_val j_nbr k_nbr c_dec_tension_chordae] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 

            tension = k_val; 

            if decreasing_tension && (k_val ~= 0.0)
                tension = tension + k_val * tension_decreasing(C(:,i), nbr, du_chordae, c_dec_tension_chordae) ; 
            end

            if tension_debug
                dec = tension_decreasing(C(:,i), nbr, du, c_dec_tension_chordae) ; 
                fprintf('tension = %e, dec_tension = %f, (i, nbr_idx, tree_idx) = (%d, %d, %d) chordae\n', tension, dec, i, nbr_idx, tree_idx); 
            end 

            % ensure that the circumferential figure is current 
            % figure(fig); 

            x_vals = [C(1,i), nbr(1)]; 
            y_vals = [C(2,i), nbr(2)]; 
            z_vals = [C(3,i), nbr(3)]; 

            % fraction of maximum tension gives fraction of way
            % through color bar 
            tension_frac_of_max = tension/max_tension;  
            color_idx = ceil(tension_frac_of_max * n_colors); 

            if color_idx == 0
                color_idx = 1; 
            end 

            if color_idx > n_colors
                %warning('max color in chordae')
                color_idx = n_colors; 
            end 

            plot3(x_vals,y_vals,z_vals,'color',cmap(color_idx,:)); 
            
            tension_by_tangent = tension * (nbr - C(:,i)) / norm(nbr - C(:,i));  

            F_chordae(tree_idx).C(:,i) = F_chordae(tree_idx).C(:,i) + tension_by_tangent; 

        end 
    end 
end 




 


