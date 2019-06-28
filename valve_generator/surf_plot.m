function fig = surf_plot(leaflet, fig)
% 
% Plots the surface and chordae
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

X_copy      = leaflet.X; 
j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc; 

if isfield(leaflet, 'name') && strcmp(leaflet.name, 'aortic')
    
    % anterior part 
    x_component = squeeze(X_copy(1,:,:)); 
    y_component = squeeze(X_copy(2,:,:)); 
    z_component = squeeze(X_copy(3,:,:)); 

    width = 1.0; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

    axis equal 
    axis auto 
    hold on 
    
else 

    % mitral default  
    j_range_anterior   = leaflet.j_range_anterior; 
    j_range_right_comm = leaflet.j_range_right_comm; 
    j_range_posterior  = leaflet.j_range_posterior; 
    j_range_left_comm  = leaflet.j_range_left_comm; 

    % NaN mask in the copy 
    for j=1:j_max
        for k=1:k_max
            if ~(is_internal(j,k) || is_bc(j,k))
               X_copy(:,j,k) = NaN;  
            end
        end 
    end

    if isfield(leaflet, 'reflect_x') && leaflet.reflect_x
        X_copy(1,:,:)     = -X_copy(1,:,:);  
    end 

    % open the figure if not passed in 
    if ~exist('fig', 'var')
        fig = figure; 
    end 

    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
        n_periodic = sum(periodic_j); 

        % anterior part 
        x_component = squeeze(X_copy(1,j_range_anterior,:)); 
        y_component = squeeze(X_copy(2,j_range_anterior,:)); 
        z_component = squeeze(X_copy(3,j_range_anterior,:)); 

        width = 1.0; 
        surf(x_component, y_component, z_component, 'LineWidth',width);

        axis equal 
        axis auto 
        hold on 

        x_component = squeeze(X_copy(1,j_range_right_comm,:)); 
        y_component = squeeze(X_copy(2,j_range_right_comm,:)); 
        z_component = squeeze(X_copy(3,j_range_right_comm,:));
        surf(x_component, y_component, z_component, 'LineWidth',width);

        x_component = squeeze(X_copy(1,j_range_posterior,:)); 
        y_component = squeeze(X_copy(2,j_range_posterior,:)); 
        z_component = squeeze(X_copy(3,j_range_posterior,:));
        surf(x_component, y_component, z_component, 'LineWidth',width);

        x_component = squeeze(X_copy(1,j_range_left_comm,:)); 
        y_component = squeeze(X_copy(2,j_range_left_comm,:)); 
        z_component = squeeze(X_copy(3,j_range_left_comm,:));
        surf(x_component, y_component, z_component, 'LineWidth',width);

        % patch commissures
        % add periodic hoops all around 
        if n_periodic > 1

            hoops = X_copy(:, [1:j_max, 1], (k_max-n_periodic+1):k_max); 

            x_component = squeeze(hoops(1,:,:)); 
            y_component = squeeze(hoops(2,:,:)); 
            z_component = squeeze(hoops(3,:,:)); 
            surf(x_component, y_component, z_component, 'LineWidth',width);

        end

    else 

        x_component = squeeze(X_copy(1,:,:)); 
        y_component = squeeze(X_copy(2,:,:)); 
        z_component = squeeze(X_copy(3,:,:)); 

        width = 1.0; 
        surf(x_component, y_component, z_component, 'LineWidth',width);

        axis equal 
        axis auto 
        hold on 

    end 



    % add chordae 
    if isfield(leaflet, 'chordae_tree') && leaflet.chordae_tree 
        for tree_idx = 1:leaflet.num_trees
            tree_plot(leaflet, tree_idx, fig); 
        end 
    end 

end 
