function fig = surf_plot(leaflet, fig)
% 
% Plots the surface and chordae
% 

X_copy      = leaflet.X; 
N           = leaflet.N; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc; 

% NaN mask in the copy 
% fill in the 3d array
% loop through N+1 to include the ring points here 

if leaflet.radial_and_circumferential
    error('plot not implemented for radial fibers') 
end

for j=1:N+1
    for k=1:N+1
        if ~(is_internal(j,k) || is_bc(j,k))
           X_copy(:,j,k) = NaN;  
        end
    end 
end

if leaflet.reflect_x
    X_copy(1,:,:)     = -X_copy(1,:,:);  
end 

% open the figure if not passed in 
if nargin < 3
    fig = figure; 
end 

x_component = squeeze(X_copy(1,:,:)); 
y_component = squeeze(X_copy(2,:,:)); 
z_component = squeeze(X_copy(3,:,:)); 

width = 1.5; 
surf(x_component, y_component, z_component, 'LineWidth',width);

axis equal 
axis auto 
hold on 

% clean up the bc on the whole surface 
string_x = [X_copy(1,1,N), X_copy(1,1,N+1)];
string_y = [X_copy(2,1,N), X_copy(2,1,N+1)];
string_z = [X_copy(3,1,N), X_copy(3,1,N+1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

string_x = [X_copy(1,N,1), X_copy(1,N+1,1)];
string_y = [X_copy(2,N,1), X_copy(2,N+1,1)];
string_z = [X_copy(3,N,1), X_copy(3,N+1,1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

% add chordae 
tree_plot(leaflet, fig); 














