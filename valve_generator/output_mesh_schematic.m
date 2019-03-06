function [] = output_mesh_schematic(valve)
    % 
    % outputs a schematic of the current mesh 
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

    pnas_figure = true; 
    
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
    alpha                  = leaflet.alpha; 
    beta                   = leaflet.beta; 

    % same number of points, crop to 2D array 
    X_schematic = nan * zeros(size(X_current)); 
    % X_schematic = nan * X_schematic(1:2,:,:); 

    plot_options = 'k-'; 
    marker_size = 3; 

    for j=1:j_max
        for k=1:k_max
            if is_internal(j, k) || is_bc(j,k)
                % X_schematic(:,j,k) = du * [j;k;nan]; 
                u = du * j;  
                if j > wrap_idx
                    % X_schematic(1,j,k) = X_schematic(1,j,k) - 1; 
                    u = u - 1;
                end
                
                v = 1 - (k_max - k) * du;
                X_schematic(:,j,k) = [u,v,nan];
            end 
        end 
    end 

    % undoing periodic wraps 
    min_x = min(min(X_schematic(1,:,:))); 
    X_schematic(1,:,:) = X_schematic(1,:,:) - min_x; 
    
    if valve.commissural_leaflets
        
        X_schematic(1,:,:) = X_schematic(1,:,:) - min(min(X_schematic(1,:,:))) + du/2; 
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
        center = .5; 
        radius = .25; 
        papillary_left_x  = center - radius + du * (-5:2:5); 
        papillary_right_x = center + radius + du * (-5:2:5);        
    end 

    papillary_x = papillary_left_x( (end - trees_anterior_left + 1) : end); 
    papillary_x = [papillary_x, papillary_right_x]; 
    papillary_x = [papillary_x, papillary_left_x(1 : (end - trees_anterior_left) )]; 
    
    % update for min_x translation
    papillary_x = papillary_x - min_x; 
    
    y_offset_trees = .45; 
    for tree_idx = 1:leaflet_schematic.num_trees
        leaflet_schematic.papillary(:,tree_idx) = [papillary_x(tree_idx); y_offset_trees; nan]; 
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

                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
 
                    if valid 
                         
                        alpha_tmp     = alpha(j_spr,k_spr); 
                        
                        x_tmp(1) = X_schematic(1,j,k); 
                        x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                        y_tmp(1) = X_schematic(2,j,k); 
                        y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 
                        
                        if alpha_tmp ~= 0
                            % lazy hack to not plot weird periodic links 
                            if abs(x_tmp(1) - x_tmp(2)) < 2*du
                                plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                            end 
                        end 

                    end 
                end 

                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 

                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                    if valid
                        beta_tmp = beta(j_spr,k_spr); 
                        x_tmp(1) = X_schematic(1,j,k); 
                        x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                        y_tmp(1) = X_schematic(2,j,k); 
                        y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 
                        if beta_tmp ~= 0
                            plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                        end
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
    if ~pnas_figure
        axis off; 
    end
    set(gcf,'color',[1 1 1])
    
    if pnas_figure
        xlabel('u')
        ylabel('v')
        
        m = axis
        m(1) = m(1) - 2*du    % xmin
        m(3) = m(3) - 2*du    % ymin
        m 
        % axis([-0.2656    0.7578   -0.0156    0.5234]); 
        axis(m)
        
        % pos_orig = get(gca, 'Position'); 
        % this approach fails, because setting 
        % set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        
%         ax = gca;          
%         pos = get(ax, 'Position');  
%         offsets = get(ax, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]; 
%         ax.ActivePositionProperty = 'position';
%         set(ax, 'OuterPosition', pos + offsets);
        
%         ax = gca;
%         outerpos = ax.OuterPosition;
%         ti = ax.TightInset; 
%         left = outerpos(1) + ti(1);
%         bottom = outerpos(2) + ti(2);
%         ax_width = outerpos(3) - ti(1) - ti(3);
%         ax_height = outerpos(4) - ti(2) - ti(4);
%         ax.Position = [left bottom ax_width ax_height];

        % fig = tightfig(fig); 

        % set(gcf,'paperpositionmode','auto')
        % print(gcf,'-depsc2','-loose','mesh_schematic_pnas.eps');
        disp('Try in r2009 if this adds whitespace.')
        printfig(fig, 'mesh_schematic_pnas')
    else
        if valve.commissural_leaflets
            printfig(fig, 'mesh_schematic_commissure')
        else 
            printfig(fig, 'mesh_schematic')
        end 
    end 
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



































