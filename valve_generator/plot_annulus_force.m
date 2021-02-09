function fig = plot_annulus_force(file_name, output_name, plot_ring, axis_range)

file = fopen(file_name, 'r'); 
data = fscanf(file, '%f'); 

stride = 1; 

dynes_to_newtons = 1e-5; 

force_scale = 1e-4 / dynes_to_newtons; 

data = reshape(data, 6, []); 

x  = data(1,1:stride:end);
y  = data(2,1:stride:end);
z  = data(3,1:stride:end);
fx = dynes_to_newtons * data(4,1:stride:end);
fy = dynes_to_newtons * data(5,1:stride:end);
fz = dynes_to_newtons * data(6,1:stride:end);

n_pts = length(x); 

fig = figure; 
hold on; 

scale = 0; 

set(fig, 'Position', [100, 100, 500, 250])
set(fig,'PaperPositionMode','auto')
set(fig, 'Renderer', 'Painters');

% oldcmap = colormap('hot');
% cmap = colormap( flipud(oldcmap) );

cmap = colormap('jet'); 

% n_colors = 500; 
% extended = true; 
% colormap(make_colormap(n_colors, extended)); 
% cmap = colormap;

% frac_to_use = 0.8; 
% max_color_idx_used = floor(frac_to_use * size(cmap,1)); 
% cmap = cmap(1:max_color_idx_used,:); 

n_colors = size(cmap,1); 
colors = zeros(n_pts,3); 

f_norms = zeros(n_pts,1); 

for i=1:n_pts
    f_norms(i) = norm([fx(i), fy(i), fz(i)]); 
end 

f_norm_max = max(f_norms); 

f_norm_frac_of_max = f_norms / f_norm_max; 

min_range = min(f_norms) %0.8e4  
max_range = max(f_norms) % 0.95e4 

for i=1:n_pts
    
    frac_of_range = (f_norms(i) - min_range) / (max_range - min_range); 
    
    color_idx = floor(frac_of_range * n_colors);
    
    if color_idx == 0
        color_idx = 1; 
    end 
    if color_idx > n_colors
        color_idx = n_colors; 
    end         
    
    colors(i,:) = cmap(color_idx,:); 
    
end 


fx = force_scale * fx;
fy = force_scale * fy;
fz = force_scale * fz; 

% quiver3(x,y,z,fx,fy,fz,scale, 'MaxHeadSize', 0.1); 
% plot each individually to get colors 
for i=1:n_pts
    quiver3(x(i), y(i), z(i), fx(i), fy(i), fz(i), scale, 'MaxHeadSize', 0.4, 'Color', colors(i,:));
end 

axis equal 
if exist('axis_range', 'var')
    % axis_default = [-1.998956798600000   1.998956798600000  -1.424871904700000   1.493493602900000  -0.619099804508580   0.240872864895460]; 
    axis(axis_range)
end
view(82,13)
axis equal 
grid on 

xticks([-1 0 1])
yticks([-1 0 1])
zticks([-0.5 0])

fontsize = 14; 
ax = gca; 
ax.FontSize = fontsize;


n_ticks = 3; 
tick_array = linspace(0,1,n_ticks); 
tick_labels = {};
tick_max = 1; 
legend_scale = 1; 
for n=1:n_ticks
    % tick=tick_array(n); 
    tmp = legend_scale * ((n-1) * (max_range - min_range)/(n_ticks-1) + min_range); 
    tick_labels{n} = sprintf('%.2e', tmp); 
end

cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels);

position_orig_axis = ax.Position; 
axis off 

% cbar.Position = [0.845160000000000 0.110000000000000 0.032000000000000 0.815000000000000]
cbar.Position = [0.8 0.25 0.032000000000000 0.4]

cbar.Label.String = {'tension',sprintf('(Newtons)', log10(1/legend_scale))};
cbar.Label.FontSize = fontsize; 
cbar.Label.Rotation = 0;
cbar.Label.Position = [0.4 1.5];

if plot_ring
    hold on 
    
    x_tmp = x; 
    y_tmp = y; 
    z_tmp = z; 
    x_tmp = [x_tmp, x_tmp(1)]; 
    y_tmp = [y_tmp, y_tmp(1)]; 
    z_tmp = [z_tmp, z_tmp(1)]; 
    plot3(x_tmp,y_tmp,z_tmp,'k')

end 

axis equal 
if exist('axis_range', 'var')
    axis(axis_range)
end 

ax.Position = position_orig_axis; 

print(fig, '-depsc', output_name);


