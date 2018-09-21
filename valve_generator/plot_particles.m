function fig = plot_particles(file_name, max_velocity, fig, bounding_box, colorbar_figure, line_width_leaflet, line_width_tails, dot_size, colored_heads)

if ~exist('fig', 'var')
    fig = figure;     
end 
% set(fig,'Renderer','painters')

if ~exist('bounding_box', 'var')
    bounding_box = false;     
end 

if ~exist('colorbar_figure', 'var')
    colorbar_figure = false;     
end

if ~exist('colorbar_for_movie', 'var')
    colorbar_for_movie = false;     
end 

if ~exist('line_width_leaflet', 'var')
    line_width_leaflet = 0.5; 
end 

if ~exist('line_width_tails', 'var')
    line_width_tails = 0.5; 
end 

if ~exist('dot_size', 'var')
    dot_size = 1.0; 
end 

if ~exist('colored_heads', 'var') 
    colored_heads = false; 
end 

run(file_name); 


plot3(x_coords_springs, y_coords_springs, z_coords_springs, '-k', 'LineWidth', line_width_leaflet); 
% plot3(x_coords_springs(1:40:end), y_coords_springs(1:40:end), z_coords_springs(1:40:end), '-k'); 
view(0,0)

hold on 




% process speeds to colormap 
all_colors = true; 
n_colors = 500;
colormap_temp = make_colormap(n_colors); 
if all_colors 
    range = 1:floor(   length(colormap_temp)); 
else 
    range = 1:floor(.9*length(colormap_temp)); 
end 
colormap_croppeed = colormap_temp(range,:); 
colormap(colormap_croppeed); 
cmap = colormap;
n_colors = size(cmap,1); 



if colorbar_figure 
    % colorbar stuff 
    n_ticks = 4; 
    tick_array = linspace(0,1,n_ticks); 
    tick_labels = {}; 
    for i=1:length(tick_array)
        tick=tick_array(i); 
        v = tick * max_velocity; 
        tick_labels{i} = sprintf('%.2f', v); 
    end 

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


if bounding_box 
    L = 3; 

    x = [-L,L]; 
    y = [-L,L]; 
    z = [3.0 - 4*L,3]; 

    % always 8 vertices on a rectangular prism 
    vertices   = zeros(3,8); 
    vertex_idx = 1; 

    for i=1:2
        for j=1:2
            for k=1:2 
                vertices(:,vertex_idx) = [x(i); y(j); z(k)]; 
                vertex_idx = vertex_idx + 1; 
            end 
        end 
    end 

    for i=1:7
        for j=(i+1):8
            current = vertices(:,i); 
            nbr = vertices(:,j); 
            % vertex is connected if at least two components are equal 
            if sum(current == nbr) >= 2
                x = [current(1), nbr(1)]; 
                y = [current(2), nbr(2)]; 
                z = [current(3), nbr(3)]; 
                plot3(x,y,z, '-k')
            end 
        end 
    end
end 


comet_tails = true; 
if comet_tails 
    velocity_frac_of_max = particle_velocity/max_velocity;
    color_idx = floor(velocity_frac_of_max * n_colors);

    % correct for maximum and minimum
    color_idx(color_idx == 0)       = 1; 
    color_idx(color_idx > n_colors) = n_colors; 

    color_idx(isnan(color_idx))     = 1; 


    % plot comet tails 
    % plot3(x_comet_coords, y_comet_coords, z_comet_coords, '-')

    comet_tail_length = size(particle_velocity,1); 
    n_particles       = size(particle_velocity,2); 

    for k = 1:n_particles
        for j = 1:(comet_tail_length-1) 
            plot3(x_comet_coords(j:j+1,k), y_comet_coords(j:j+1,k), z_comet_coords(j:j+1,k), '-', 'color', cmap(color_idx(j,k),:), 'LineWidth', line_width_tails); 
        end 
    end 



    cheap_comet_heads = true; 
    if cheap_comet_heads
        % comet heads -- cheap 

        if colored_heads
            scatter3(x_comet_coords(1,:), y_comet_coords(1,:), z_comet_coords(1,:), dot_size, cmap(color_idx(1,:),:), 'filled') % weird behavior on 'filled option here'
        else
            scatter3(x_comet_coords(1,:), y_comet_coords(1,:), z_comet_coords(1,:), dot_size, 'k', 'filled') % weird behavior on 'filled option here'
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

end 

axis equal; 
axis tight; 
axis off; 
