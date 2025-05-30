function [coapt_height, coapt_reserve_height, belly_height, coapt_min, coapt_max, fig, min_distances] = coaptation_analysis_aortic(leaflets, threshold, plots, fig)

j_max  = leaflets(1,1).j_max; 
k_max  = leaflets(1,1).k_max; 
N_each = leaflets(1,1).N_each; 

[n_layers, n_leaflets] = size(leaflets);


is_internal = leaflets(1,1).is_internal; 
is_bc = leaflets(1,1).is_bc; 
% N_each = leaflets(1,1).N_each; 

layer_idx_ventricular = 1; 

if ~exist('fig', 'var')
    if plots 
        fig = figure;  
    else 
        fig = [];
    end 
end 

debug_plots = false; 
if debug_plots
    load plot_data_temp.mat
else 

    if isfield(leaflets(1,1), 'N_leaflets')
        if leaflets(1,1).N_leaflets ~= n_leaflets
            error('inconsistent number of leaflets')
        end         
    end 

    if n_leaflets ~= 2
        warning('not implmented, two leaflets required');     
    end
    
    min_distances = nan(n_leaflets,j_max,k_max);

    for leaflet_idx = 1:n_leaflets 
        
        if leaflet_idx == 1
            leaflet_comparison_idx = 2;
        else 
            leaflet_comparison_idx = 1;
        end 
        
        % twod arrays 
        x_comp_comparison = leaflets(layer_idx_ventricular, leaflet_comparison_idx).X(1,2:(j_max-1), :); 
        y_comp_comparison = leaflets(layer_idx_ventricular, leaflet_comparison_idx).X(2,2:(j_max-1), :); 
        z_comp_comparison = leaflets(layer_idx_ventricular, leaflet_comparison_idx).X(3,2:(j_max-1), :); 

        % flattened 
        x_comp_comparison = x_comp_comparison(:); 
        y_comp_comparison = y_comp_comparison(:); 
        z_comp_comparison = z_comp_comparison(:); 

        triangulation = delaunay(x_comp_comparison, y_comp_comparison, z_comp_comparison); 

        X_comparison_matrix = [x_comp_comparison, y_comp_comparison, z_comp_comparison]; 



        for j=1:j_max
        %    for k=1:k_max         
            X_row_temp = squeeze(leaflets(layer_idx_ventricular, leaflet_idx).X(:,j,:))';
            [~, min_distances(leaflet_idx,j,:)] = dsearchn(X_comparison_matrix, triangulation, X_row_temp);                 
        %    end     
        end 

    end 
    
    coapted_points = (min_distances < threshold); 

end 



% coaptation height and reserve height  
coapt_height = 0.0;
coapt_reserve_height = 0.0; 
coapt_reserve_initialized = false; 
belly_height = 0.0; 



leaflet_idx = 1; 

center_idx_j = N_each/2; 
j = center_idx_j; 

X_min_coapted = nan(3,1); 

for k=2:k_max

    j_nbr_tmp = j;
    k_nbr_tmp = k-1; 
    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflets(layer_idx_ventricular, leaflet_idx), j, k, j_nbr_tmp, k_nbr_tmp); 
    if ~valid 
        error('trying to compute lengths with an invalid rest length')
    end

    X_temp = leaflets(layer_idx_ventricular, leaflet_idx).X(:,j,k);
    X_nbr = leaflets(layer_idx_ventricular, leaflet_idx).X(:,j_nbr,k_nbr); 
    
    % first relevant point 
    if coapted_points(leaflet_idx,j,k) && (~coapt_reserve_initialized)
        X_min_coapted = X_temp;
    end 

    if coapted_points(leaflet_idx,j,k) && coapted_points(leaflet_idx,j,k-1)
        coapt_height = coapt_height + norm(X_temp - X_nbr);
        coapt_reserve_initialized = true; 
    end 

    if coapt_reserve_initialized
        coapt_reserve_height = coapt_reserve_height + norm(X_temp - X_nbr);
    else 
        belly_height = belly_height + norm(X_temp - X_nbr);
    end 
    
end 
    
X_min_coapted;
X_max_coapted = leaflets(layer_idx_ventricular, leaflet_idx).X(:,center_idx_j,k_max);

coapt_min = X_min_coapted(3); 
coapt_max = leaflets(layer_idx_ventricular, leaflet_idx).X(3,center_idx_j,k_max);     


if plots
       
    max_plot_cap = max(max(max(min_distances))); 

    n_colors = 500;
    extended = false; 
    colormap(flipud(make_colormap(n_colors, extended))); 
    cmap = colormap;
    n_colors = size(cmap,1); 

    % plot the actual surface 
    X_copy      = leaflets(layer_idx_ventricular, leaflet_idx).X; 

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

    outline_cleanup = true; 
    if outline_cleanup 

        % horizontal at bottom and free edge 
        for k=[1,k_max]
            for j=1:N_each

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
            for j=[1,j_max]

                X = X_copy(:,j,k); 
                X_nbr = X_copy(:,j,k+1); 

                x_vals = [X(1), X_nbr(1)]; 
                y_vals = [X(2), X_nbr(2)]; 
                z_vals = [X(3), X_nbr(3)]; 

                plot3(x_vals,y_vals,z_vals,'k'); 
            end 
        end 
    end 


    tick_max = max_plot_cap; 
    colors = zeros(j_max,k_max,3); 
    for j=1:j_max
        for k=1:k_max
            color_idx = floor(n_colors * min_distances(leaflet_idx,j,k) / max_plot_cap); 
            if (color_idx == 0) || (isnan(color_idx))
                color_idx = 1; 
            end 
            if color_idx > n_colors
                color_idx = n_colors; 
            end         
            colors(j,k,:) = cmap(color_idx,:); 
        end 
    end 

%     circ_off_by_one_adjust = true; 
%     if (circ || ratio) && circ_off_by_one_adjust 
%         % colors are read from the minimum j,k index for the square 
%         % so shift down by one 
%         colors_temp = zeros(size(colors)); 
% 
%         for j=1:j_max 
%             for k=1:(k_max-1)
%                 colors_temp(j,k,:) = colors(j,k+1,:); 
%             end 
%         end 
% 
%         colors = colors_temp; 
%     end 

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
        tick_labels{n} = sprintf('%.2f', tmp); 
    end 

    colorbar_on = true; 
    if colorbar_on
        colorbar('Ticks', tick_array, 'TickLabels', tick_labels);
    end 

    axis equal 
    axis off 
    axis tight

end 












