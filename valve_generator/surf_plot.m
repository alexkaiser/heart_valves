function fig = surf_plot(leaflet, fig)
% 
% Plots the surface and chordae
% 

X_copy      = leaflet.X; 
N           = leaflet.N; 
j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc; 

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
    
    N_anterior = j_max/2; 
    
    % anterior part 
    x_component = squeeze(X_copy(1,1:N_anterior,:)); 
    y_component = squeeze(X_copy(2,1:N_anterior,:)); 
    z_component = squeeze(X_copy(3,1:N_anterior,:)); 

    width = 1.5; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

    axis equal 
    axis auto 
    hold on 

    
    x_component = squeeze(X_copy(1,(N_anterior+1):j_max,:)); 
    y_component = squeeze(X_copy(2,(N_anterior+1):j_max,:)); 
    z_component = squeeze(X_copy(3,(N_anterior+1):j_max,:));
    surf(x_component, y_component, z_component, 'LineWidth',width);
    
    % patch commissures
    % ring is always periodic and do not patch it 
    if n_periodic > 1

        left_comm = X_copy(:, [j_max, 1], (k_max-n_periodic+1):k_max); 

        x_component = squeeze(left_comm(1,:,:)); 
        y_component = squeeze(left_comm(2,:,:)); 
        z_component = squeeze(left_comm(3,:,:)); 
        surf(x_component, y_component, z_component, 'LineWidth',width);

        right_comm = X_copy(:,N_anterior:(N_anterior+1),(k_max-n_periodic+1):k_max); 

        x_component = squeeze(right_comm(1,:,:)); 
        y_component = squeeze(right_comm(2,:,:)); 
        z_component = squeeze(right_comm(3,:,:)); 
        surf(x_component, y_component, z_component, 'LineWidth',width); 
            
    end
        
else 
    
    x_component = squeeze(X_copy(1,:,:)); 
    y_component = squeeze(X_copy(2,:,:)); 
    z_component = squeeze(X_copy(3,:,:)); 

    width = 1.5; 
    surf(x_component, y_component, z_component, 'LineWidth',width);

    axis equal 
    axis auto 
    hold on 

end 


% if isfield(leaflet, 'radial_and_circumferential')
%     if leaflet.radial_and_circumferential
% 
%         % clean up the bc on the whole surface 
%         string_x = [X_copy(1,1,k_max), X_copy(1,2,k_max)];
%         string_y = [X_copy(2,1,k_max), X_copy(2,2,k_max)];
%         string_z = [X_copy(3,1,k_max), X_copy(3,2,k_max)];
%         plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
% 
%         string_x = [X_copy(1,j_max-1,k_max), X_copy(1,j_max,k_max)];
%         string_y = [X_copy(2,j_max-1,k_max), X_copy(2,j_max,k_max)];
%         string_z = [X_copy(3,j_max-1,k_max), X_copy(3,j_max,k_max)];
%         plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
% 
%     else
% 
%         % clean up the bc on the whole surface 
%         string_x = [X_copy(1,1,N), X_copy(1,1,N+1)];
%         string_y = [X_copy(2,1,N), X_copy(2,1,N+1)];
%         string_z = [X_copy(3,1,N), X_copy(3,1,N+1)];
%         plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
% 
%         string_x = [X_copy(1,N,1), X_copy(1,N+1,1)];
%         string_y = [X_copy(2,N,1), X_copy(2,N+1,1)];
%         string_z = [X_copy(3,N,1), X_copy(3,N+1,1)];
%         plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
% 
%     end 
% end 

% add chordae 
if isfield(leaflet, 'chordae_tree') && leaflet.chordae_tree 
    for tree_idx = 1:leaflet.num_trees
        tree_plot(leaflet, tree_idx, fig); 
    end 
end 


