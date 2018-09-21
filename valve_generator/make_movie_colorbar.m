function make_movie_colorbar(max_velocity, format_string, time, file_name)

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


% colorbar for movie, bigger and bolder 
n_ticks = 4; 
tick_array = linspace(0,1,n_ticks); 
tick_labels = {}; 
for i=1:length(tick_array)
    tick=tick_array(i); 
    v = tick * max_velocity; 
    tick_labels{i} = sprintf('%d', v); 
end 

fontsize = 24 * 4;

fig = figure('visible','off'); 
set(fig, 'Renderer', 'Painters');
set(fig, 'Position', [0 0 444 4320])
set(fig, 'PaperPositionMode','auto')

colormap(colormap_croppeed); 
cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels, 'fontweight', 'bold', 'fontsize',fontsize); 
cbar.Label.String = sprintf('time (s)\n%.2f\n\n\n|u|\n(cm/s)', time);
cbar.Label.FontSize = fontsize; 
cbar.Label.Position = [0 1.05]; 
cbar.Label.Rotation = 0;
cbar.Label.FontWeight = 'bold'; 
cbar.TickDirection = 'out'; 
cbar.AxisLocation = 'in'; 


% position and sizing of the bar itself relative to the block 
% does not include labels
cbar.Position = [0.6 0.27 .4 .5]; 


grid off 
axis off 

print(fig, format_string, file_name, '-r0'); 


