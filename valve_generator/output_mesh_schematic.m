function [] = output_mesh_schematic(valve)
    % 
    % outputs a schematic of the current mesh 
    % 

    leaflet = valve.leaflets(1); 

    % move points in periodic way past this index
    % hacks to make mesh shape be split in loction of dissection 
    wrap_idx = valve.N_anterior + valve.N_posterior/2; 
    if valve.commissural_leaflets
        %error('not implemented for schematic, must manually for papillary arrangement'); 
        wrap_idx = wrap_idx + valve.N_commissure; 
    end 

    X_current              = leaflet.X; 
    j_max                  = leaflet.j_max; 
    k_max                  = leaflet.k_max; 
    du                     = leaflet.du; 
    is_internal            = leaflet.is_internal; 
    is_bc                  = leaflet.is_bc; 


    % same number of points, crop to 2D array 
    X_schematic = nan * zeros(size(X_current)); 
    % X_schematic = nan * X_schematic(1:2,:,:); 

    plot_options = 'k-'; 
    marker_size = 3; 

    for j=1:j_max
        for k=1:k_max
            if is_internal(j, k) || is_bc(j,k)
                X_schematic(:,j,k) = du * [j;k;nan]; 
                
                if j > wrap_idx
                    X_schematic(1,j,k) = X_schematic(1,j,k) - 1; 
                end 
                
            end 
        end 
    end 

    leaflet_schematic   = leaflet; 
    leaflet_schematic.X = X_schematic; 
        

    % crop for 2d papillary 
    % leaflet_schematic.papillary = leaflet_schematic.papillary(1:2,:); 

    % left and right papillary hack placement 
    if ~valve.commissural_leaflets
        trees_anterior_left = 1; 
        % trees_per_side = size(leaflet.papillary,2)/2; 
        papillary_left_x  =      du * [-3, -2, -1, 1]; 
        papillary_right_x = .5 + du * [-1,  1,  2, 3]; 
    else 
        trees_anterior_left = 2; 
        % trees_per_side = size(leaflet.papillary,2)/2; 
        papillary_left_x  =      du * (-5:0); 
        papillary_right_x = .5 + du * ( 0:5);        
    end 

    papillary_x = papillary_left_x( (end - trees_anterior_left + 1) : end); 
    papillary_x = [papillary_x, papillary_right_x]; 
    papillary_x = [papillary_x, papillary_left_x(1 : (end - trees_anterior_left) )]; 
    
    for tree_idx = 1:leaflet_schematic.num_trees
        leaflet_schematic.papillary(:,tree_idx) = [papillary_x(tree_idx); 0; nan]; 
    end 
    
    for tree_idx = 1:leaflet_schematic.num_trees
        leaflet_schematic = add_chordae(leaflet_schematic, tree_idx); 
    end 


    fig = figure; 
    hold on; 

    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) 

                % plot all connections to current location 
                % redundantly becuase who cares 

                % u type fibers 
                for j_nbr_tmp = [j-1,j+1]

                    k_nbr_tmp = k; 

                    [valid j_nbr k_nbr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
 
                    if valid 
                        x_tmp(1) = X_schematic(1,j,k); 
                        x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                        y_tmp(1) = X_schematic(2,j,k); 
                        y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 
                        
                        % lazy hack to not plot weird periodic links 
                        if abs(x_tmp(1) - x_tmp(2)) < 2*du
                            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                        end 

                    end 
                end 

                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 

                    [valid j_nbr k_nbr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                    if valid
                        x_tmp(1) = X_schematic(1,j,k); 
                        x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                        y_tmp(1) = X_schematic(2,j,k); 
                        y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 
                        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                    end 
                end 
            end       
        end 
    end 
    
    % more hacks, clean up periodic boundary 
    j=wrap_idx; 
    for k=1:k_max
        if is_internal(j,k) 
            x_tmp(1) = X_schematic(1,j,k); 
            x_tmp(2) = X_schematic(1,j,k) + du/2; 
            y_tmp(1) = X_schematic(2,j,k); 
            y_tmp(2) = X_schematic(2,j,k); 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
            % make a little zig zag to show periodicity 
            x_tmp(1) = X_schematic(1,j,k) + du/2; 
            x_tmp(2) = X_schematic(1,j,k) + 3*du/4; 
            y_tmp(1) = X_schematic(2,j,k); 
            y_tmp(2) = X_schematic(2,j,k) + du/4; 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
            x_tmp(1) = X_schematic(1,j,k) + 3*du/4; 
            x_tmp(2) = X_schematic(1,j,k) + du; 
            y_tmp(1) = X_schematic(2,j,k) + du/4; 
            y_tmp(2) = X_schematic(2,j,k); 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
        end 
    end 
    
    j=wrap_idx + 1; 
    for k=1:k_max
        if is_internal(j,k) 
            x_tmp(1) = X_schematic(1,j,k); 
            x_tmp(2) = X_schematic(1,j,k) - du/2; 
            y_tmp(1) = X_schematic(2,j,k); 
            y_tmp(2) = X_schematic(2,j,k); 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
            % make a little zig zag to show periodicity 
            x_tmp(1) = X_schematic(1,j,k) - du/2; 
            x_tmp(2) = X_schematic(1,j,k) - 3*du/4; 
            y_tmp(1) = X_schematic(2,j,k); 
            y_tmp(2) = X_schematic(2,j,k) - du/4; 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
            x_tmp(1) = X_schematic(1,j,k) - 3*du/4; 
            x_tmp(2) = X_schematic(1,j,k) - du; 
            y_tmp(1) = X_schematic(2,j,k) - du/4; 
            y_tmp(2) = X_schematic(2,j,k); 
            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
            
        end 
    end 
    
    
    
    for tree_idx = 1:leaflet.num_trees
        tree_plot_flat(leaflet_schematic, tree_idx, fig, plot_options, marker_size); 
    end 

    % axis([-.3, .8, -.1, .6])
    axis tight; 
    axis equal;    
    axis off; 
    set(gcf,'color',[1 1 1])
    
    % printfig(fig, 'mesh_schematic')
    
end 



function [] = tree_plot_flat(leaflet, tree_idx, fig, plot_options, marker_size)
    % 
    % Plots chordae tendineae on figure 'fig'
    % 
    % leaflet     Big data strucure 
    % fig         Current figure handle 

    set(0, 'currentfigure', fig); 
    hold on 

    % unpack some data 
    X              = leaflet.X; 
    root           = leaflet.chordae(tree_idx).root; 
    free_edge_idx  = leaflet.chordae(tree_idx).free_edge_idx; 
    C              = leaflet.chordae(tree_idx).C; 
 

    [n_leaves, m] = size(free_edge_idx);  
    n_tree = log2(n_leaves);

    % there are max_internal points on the tree 
    % leaves are not included as they are part of the leaflet 
    [m max_internal] = size(C); 

    if m ~= 3
        error('chordae must be three dimensional')
    end 

    % this parameter is the same N as in the leaflet 
    % it is the number of (not included) leaves in the tree
    if n_leaves ~= (max_internal + 1)
        error('inconsistent dimensions of chordae and free edge indices')
    end 

    % sanity check in building a balanced tree 
    if abs(n_tree - round(n_tree)) > eps 
        error('must use a power of two'); 
    end 

    % left side 
    x = [root(1); C(1,1)]; 
    y = [root(2); C(2,1)]; 
    plot(x,y,plot_options,'MarkerSize', marker_size); 


    for i = 1:max_internal

        left_child_idx = 2*i; 

        if left_child_idx <= max_internal
            left_child = C(:,left_child_idx); 
        else 
            j = left_child_idx - max_internal; 
            left_child = X(:, free_edge_idx(j,1), free_edge_idx(j,2));
        end 

        x = [C(1,i); left_child(1)]; 
        y = [C(2,i); left_child(2)]; 
        plot(x,y,plot_options, 'MarkerSize', marker_size); 

        right_child_idx = 2*i + 1; 

        if right_child_idx <= max_internal
            right_child = C(:,right_child_idx); 
        else 
            j = right_child_idx - max_internal; 
            right_child = X(:, free_edge_idx(j,1), free_edge_idx(j,2));
        end 

        x = [C(1,i); right_child(1)]; 
        y = [C(2,i); right_child(2)]; 
        plot(x,y,plot_options, 'MarkerSize', marker_size); 

    end 

end 



































