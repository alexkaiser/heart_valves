function fig = surf_plot_commissure(params, filter_params, left, fig)
% 
% Plots the surface and connective chordae
% 
% 

X_copy = params.X; 
N      = params.N; 

if left 
    papillary = [0; -filter_params.a; 0]; 
else 
    papillary = [0;  filter_params.a; 0]; 
end 

if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end 

% NaN mask in the copy 

% fill in the 3d array
% loop through N+1 to include the ring points here 
for j=1:N+2
    for k=1:((N+3)/2)
        % out of the domain? 
        % apply a NaN mask for plotting 
        if ~in_domain_with_bc(j,k,N)
            X_copy(:,j,k) = NaN; 
        end

    end 
end

% open the figure if not passed in 
if nargin < 4
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

% clean up the boundaries which are ignored by surf 
string_x = [X_copy(1,1,1), X_copy(1,2,1)];
string_y = [X_copy(2,1,1), X_copy(2,2,1)];
string_z = [X_copy(3,1,1), X_copy(3,2,1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

string_x = [X_copy(1,N+1,1), X_copy(1,N+2,1)];
string_y = [X_copy(2,N+1,1), X_copy(2,N+2,1)];
string_z = [X_copy(3,N+1,1), X_copy(3,N+2,1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

string_x = [X_copy(1,(N+3)/2,(N+3)/2 - 1), X_copy(1,(N+3)/2,(N+3)/2)];
string_y = [X_copy(2,(N+3)/2,(N+3)/2 - 1), X_copy(2,(N+3)/2,(N+3)/2)];
string_z = [X_copy(3,(N+3)/2,(N+3)/2 - 1), X_copy(3,(N+3)/2,(N+3)/2)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 


% add chordae as line segments 
k = 1; 
for j=1:(N+2)
    if is_internal_commissure(j,k,N)
        string_x = [papillary(1), X_copy(1,j,k)];
        string_y = [papillary(2), X_copy(2,j,k)];
        string_z = [papillary(3), X_copy(3,j,k)];
        plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
    end 
end 















