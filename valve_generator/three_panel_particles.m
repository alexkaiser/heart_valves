

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

    extra = .2; 
    L = 3; 
    horiz_min = -L - extra; 
    horiz_max =  L + extra; 
    zmin =  3 - 4*L - extra;
    zmax =  3 + extra; 
    
    file_name = strcat(data_dir, '/', file_name)

    fig_subplots = figure; 
    % 4x HD resolution 
    % set(fig_subplots, 'Position', [100 100 7680 4320])
    
    hSub1 = subplot(1,3,1);

%     hSub2 = subplot(1,3,2);
%     hSub3 = subplot(1,3,3);

    fig_subplots = plot_particles(file_name, max_velocity, fig_subplots, colorbar_figure); 
    
    % side view, anterior leaflet to the right 
    view(0,0)

    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 
    width  = horiz_max - horiz_min; 
    height = zmax - zmin; 
    
    % front view 
    
    hSub2 = subplot(1,3,2);
    fig_subplots = plot_particles(file_name, max_velocity, fig_subplots, colorbar_figure); 
    view(90,0)
    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])
    axis off 
    
    
    
    
    
    
    
    
    
    
    % unclear what this does, but consider if there are problems 
    % set(fig,'PaperPositionMode','auto')
    
    print(fig_subplots, '-djpeg', strcat(output_base_name, '_three_panel'))

end 
