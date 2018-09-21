

data_dir = './' 
% data_dir =  '/Users/alex/mitral_fully_discrete/plots_and_session_files/mitral_1647284_512_git_5855699_two_leaflet_mean_reverting'; 

names = ["mitral_tree_512_001005.m";]

output_names = ["three_particle_views_001005";]; 

max_velocity = 150; 


min_frame = 0; 
max_frame = 1440; 

colorbar_name = 'colorbar_movie.jpg'; 

PATH = getenv('PATH');
setenv('PATH', [PATH ':/usr/local/bin']);

eps_format = false; 
if eps_format 
    format_extension = '.eps'; 
    format_string = '-depsc'; 
else 
    format_extension = '.jpg'; 
    format_string = '-djpeg'; 
end 

make_movie_colorbar(max_velocity, format_extension, format_string); 

if min_frame ~= 0
    error('assumes minimum frame is one'); 
end

line_width_leaflet = 2; 
line_width_tails = 16; 
dot_size = 64; 
colored_heads = true; 

cleanup = true; 

rng('shuffle');


    
range = randperm(max_frame + 1) - 1;

for i=range

    file_name = sprintf('mitral_tree_512_%.6d.m', i); 
    output_base_name = sprintf('particle_views_%.6d', i); 

    start_name = sprintf('%s_start.txt', output_base_name); 
    end_name = sprintf('%s_complete.txt', output_base_name); 

    fprintf('  Checking file_name %s...\n', file_name)

    if exist(start_name, 'file') || exist(end_name, 'file') 
        fprintf('Start or end file exists, moving on\n')
        continue; 
    end      

%         if cleanup && exist(end_name, 'file')
%             fprintf('On cleanup, end file exists, moving on\n')
%             continue; 
%         end 

    system(sprintf('touch %s', start_name)); 

    fprintf('Writing frames assoc with file_name %s\n', file_name)

    clear_fucntion_cache = false; 

    if clear_fucntion_cache 

        working_dir = pwd; 
        cd(data_dir); 

        file_name_copy = file_name; 
        clear file_name
        file_name = file_name_copy;

        cd(working_dir)
    end 

    bounding_box = true; 
    colorbar_figure = false; 
    colorbar_for_movie = true; 

    L = 3; 
    horiz_min = -L; 
    horiz_max =  L; 
    zmin =  3 - 4*L;
    zmax =  3; 

    file_name = strcat(data_dir, '/', file_name); 

    fig = figure('visible','off'); 
    set(fig, 'Renderer', 'painters');
    % 4x HD resolution 
    set(fig, 'Position', [0 0 2160 4320])
    set(fig, 'PaperPositionMode','auto')

    fig = plot_particles(file_name, max_velocity, fig, bounding_box, colorbar_figure, line_width_leaflet, line_width_tails, dot_size, colored_heads); 

    % side view, anterior leaflet to the right 
    view(0,0)

    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 

    ax = gca;
    ax.Position = [0.01 0.01 .98 .98]; 

    side_name = strcat(output_base_name, '_side', format_extension); 

    print(fig, format_string, side_name, '-r0')

    % front view 

    view(90,0)
    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 

    ax = gca;
    ax.Position = [0.01 0.01 .98 .98]; 

    front_name = strcat(output_base_name, '_front', format_extension); 
    print(fig, format_string, front_name, '-r0')

    % top view
    r = 0.1 * 21.885241; 
    view(0,90)

    frac_to_right = .35; 

    axis([-r frac_to_right*r -r r -4 2])

    % set(fig, 'Position', [100, 100, 500, floor(500*(1 + frac_to_right)/2)])
    set(fig, 'Position', [0 0 2916 4320])

    ax = gca;
    ax.Position = [0.00 0.00 1 1]; 

    top_name = strcat(output_base_name, '_top', format_extension); 
    print(fig, format_string, strcat(output_base_name, '_top'), '-r0')

    close(fig); 

    % append the existing panels 
    all_name = strcat(output_base_name, '_panels', format_extension); 
    command = sprintf('convert %s %s %s %s +append %s', colorbar_name, side_name, front_name, top_name, all_name);  
    system(command,'-echo') 

    if eps_format
        all_name_eps = strcat(output_base_name, '_panels', '.eps');
        command = sprintf('convert %s %s', all_name, all_name_eps);  
        system(command,'-echo')
    end 

    system(sprintf('touch %s', end_name)); 
    system(sprintf('rm %s', start_name)); 

end 

