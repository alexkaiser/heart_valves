function fig = surf_plot(params, filter_params, fig)
% 
% Plots the surface and connective chordae
% 
% 

X_copy = params.X; 
N      = params.N; 

left_papillary = [0;-filter_params.a;0]; 
right_papillary = [0; filter_params.a;0]; 

% NaN mask in the copy 

% fill in the 3d array
% loop through N+1 to include the ring points here 
for j=1:N+1
    for k=1:N+1

        % out of the triangle? 
        % apply a NaN mask for plotting 
        if ((j+k) > (N+2))
            X_copy(:,j,k) = NaN; 
        end

    end 
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



% add chordae as line segments 
j = 1; 
for k=1:N
    string_x = [left_papillary(1), X_copy(1,j,k)];
    string_y = [left_papillary(2), X_copy(2,j,k)];
    string_z = [left_papillary(3), X_copy(3,j,k)];
    plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
end 

k = 1; 
for j=1:N
    string_x = [right_papillary(1), X_copy(1,j,k)];
    string_y = [right_papillary(2), X_copy(2,j,k)];
    string_z = [right_papillary(3), X_copy(3,j,k)];
    plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
end 















