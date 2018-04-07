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

% plot comet tails 
plot3(x_comet_coords, y_comet_coords, z_comet_coords, '-r')

cheap_comet_heads = true; 
if cheap_comet_heads
    % comet heads -- cheap 
    scatter3(x_comet_coords(1,:), y_comet_coords(1,:), z_comet_coords(1,:), 1, 'k', 'filled') % weird behavior on 'filled option here'
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
