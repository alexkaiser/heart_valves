function fig = torus_plot(torus, fig)
% 
% Plots the surface and chordae
% 

X_copy      = torus.X; 

% open the figure if not passed in 
if ~exist('fig', 'var')
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

