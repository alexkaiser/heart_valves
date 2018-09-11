

data_dir = '/Users/alex/mitral_fully_discrete/plots_and_session_files/mitral_1647284_512_git_5855699_two_leaflet_mean_reverting'; 

names = ["mitral_tree_512_001005.m";]

output_names = ["three_particle_views_001005";]; 



for i=1:length(names)
    file_name = names(i)
    output_base_name = output_names(i)

    max_velocity = 150; 

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

    
    L = 3; 
    horiz_min = -L; 
    horiz_max =  L; 
    zmin =  3 - 4*L;
    zmax =  3; 
    
    file_name = strcat(data_dir, '/', file_name)

    fig = figure; % figure('visible','off'); 
    set(fig, 'Renderer', 'Painters');
    % 4x HD resolution 
    set(fig, 'Position', [0 0 2160 4320])
    set(fig, 'PaperPositionMode','auto')
    
    fig = plot_particles(file_name, max_velocity, fig, bounding_box, colorbar_figure); 
    
    % side view, anterior leaflet to the right 
    view(0,0)

    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 

    ax = gca;
    ax.Position = [0.05 0.05 .9 .9]; 
    
    print(fig, '-djpeg', strcat(output_base_name, '_side'), '-r0')
    
    % front view 
    
    view(90,0)
    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 
        
    ax = gca;
    ax.Position = [0.05 0.05 .9 .9]; 
    
    print(fig, '-djpeg', strcat(output_base_name, '_front'), '-r0')

    % top view
    r = 0.1 * 21.885241; 
    view(0,90)

    frac_to_right = .35; 

    axis([-r frac_to_right*r -r r -4 2])

    % set(fig, 'Position', [100, 100, 500, floor(500*(1 + frac_to_right)/2)])
    set(fig, 'Position', [0 0 2916 4320])
    
    ax = gca;
    ax.Position = [0.05 0.05 .9 .9]; 
    
    print(fig, '-djpeg', strcat(output_base_name, '_top'), '-r0')
    
end 
