function fig = plot_annulus_force(file_name, plot_ring, axis_range)

file = fopen(file_name, 'r'); 
data = fscanf(file, '%f'); 

stride = 1; 

force_scale = 1e-4; 

data = reshape(data, 6, []); 

x  = data(1,1:stride:end);
y  = data(2,1:stride:end);
z  = data(3,1:stride:end);
fx = force_scale * data(4,1:stride:end);
fy = force_scale * data(5,1:stride:end);
fz = force_scale * data(6,1:stride:end);

fig = figure; 

scale = 0; 

set(fig, 'Position', [100, 100, 500, 250])
set(fig,'PaperPositionMode','auto')
set(fig, 'Renderer', 'Painters');

quiver3(x,y,z,fx,fy,fz,scale, 'MaxHeadSize', 0.1); 

axis equal 
if exist('axis_range', 'var')
    axis(axis_range)
end
view(28,26)
axis equal 

if plot_ring
    hold on 
    
    x_tmp = x; 
    y_tmp = y; 
    z_tmp = z; 
    x_tmp = [x_tmp, x_tmp(1)]; 
    y_tmp = [y_tmp, y_tmp(1)]; 
    z_tmp = [z_tmp, z_tmp(1)]; 
    plot3(x_tmp,y_tmp,z_tmp,'k')

    axis equal 
    if exist('axis_range', 'var')
        axis(axis_range)
    end 
end 
