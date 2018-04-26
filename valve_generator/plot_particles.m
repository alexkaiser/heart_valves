function fig = plot_particles(file_name, max_velocity, fig)

if ~exist('fig', 'var')
    fig = figure;     
end 
set(fig,'Renderer','painters')

run(file_name); 

plot3(x_coords_springs, y_coords_springs, z_coords_springs, '-k')
view(0,0)
axis equal; 
axis tight; 
axis off; 

hold on 


% process speeds to colormap 
n_colors = 500;
colormap_temp = make_colormap(n_colors); 
range = 1:floor(.9*length(colormap_temp)); 
colormap_croppeed = colormap_temp(range,:); 
colormap(colormap_croppeed); 
cmap = colormap;
n_colors = size(cmap,1); 

% colorbar stuff 
n_ticks = 4; 
tick_array = linspace(0,1,n_ticks); 
tick_labels = {}; 
for i=1:length(tick_array)
    tick=tick_array(i); 
    v = tick * max_velocity; 
    tick_labels{i} = sprintf('%.2f', v); 
end 
colorbar_figure = true; 
if colorbar_figure 
    fig_colorbar = figure; 
    
    colormap(colormap_croppeed); 
    cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels); 
    cbar.Label.String = 'Velocity magnitude (cm/s)';
    grid off 
    axis off 
    printfig(fig_colorbar, 'colorbar_only_comet_tails'); 
    close(fig_colorbar); 
    
    % reset current figure 
    figure(fig); 
end



velocity_frac_of_max = particle_velocity/max_velocity;
color_idx = floor(velocity_frac_of_max * n_colors);

color_idx(color_idx == 0)       = 1; 
color_idx(color_idx > n_colors) = n_colors; 

color_idx(isnan(color_idx))     = 1; 


% plot comet tails 
% plot3(x_comet_coords, y_comet_coords, z_comet_coords, '-')

comet_tail_length = size(particle_velocity,1); 
n_particles       = size(particle_velocity,2); 

for k = 1:n_particles
    for j = 1:(comet_tail_length-1) 
        plot3(x_comet_coords(j:j+1,k), y_comet_coords(j:j+1,k), z_comet_coords(j:j+1,k), '-', 'color', cmap(color_idx(j,k),:)); 
    end 
end 



cheap_comet_heads = true; 
if cheap_comet_heads
    % comet heads -- cheap 
    
    colored_heads = false; 
    if colored_heads
        scatter3(x_comet_coords(1,:), y_comet_coords(1,:), z_comet_coords(1,:), 1, cmap(color_idx(1,:),:), 'filled') % weird behavior on 'filled option here'
    else
        scatter3(x_comet_coords(1,:), y_comet_coords(1,:), z_comet_coords(1,:), 1, 'k', 'filled') % weird behavior on 'filled option here'
    end 
else 
    % this is realllllllly slow and makes 
    % source -- https://www.mathworks.com/matlabcentral/answers/254961-3d-plot-points-as-spheres-instead-of-dots
    x_heads = x_comet_coords(1,:); 
    y_heads = y_comet_coords(1,:); 
    z_heads = z_comet_coords(1,:); 

    R = .02; 
    NumSphFaces = 15;
    [SX,SY,SZ] = sphere(NumSphFaces);
    for K = 1:length(x_heads)
      surf(SX*R + x_heads(K), SY*R + y_heads(K), SZ*R + z_heads(K));
    end
    lighting phong

end 



