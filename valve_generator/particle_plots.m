

data_dir = '/Users/alex/mitral_fully_discrete/plots_and_session_files/mitral_1647284_512_git_5855699_two_leaflet_mean_reverting'; 

names = ["mitral_tree_512_001005.m"; "mitral_tree_512_001267.m"; "mitral_tree_512_001382.m"]

output_names = ["particle_views_001005";"particle_views_001267";"particle_views_001382"]; 



for i=1:length(names)
    file_name = names(i)
    output_base_name = output_names(i)

    max_velocity = 150; 

    clear_fucntion_cache = true; 

    if clear_fucntion_cache 

        working_dir = pwd; 
        cd(data_dir); 

        file_name_copy = file_name; 
        clear file_name
        file_name = file_name_copy;

        cd(working_dir)
    end 


    file_name = strcat(data_dir, '/', file_name)


    fig = figure; 
    fig = plot_particles(file_name, max_velocity); 
    % printfig(fig, strcat(output_base_name, '_side_full'))



    view(0,0)

    horiz_min = -3; 
    horiz_max =  3; 
    zmin = -7;
    zmax =  2; 

    axis([horiz_min horiz_max horiz_min horiz_max zmin zmax])


    width  = horiz_max - horiz_min; 
    height = zmax - zmin; 



    set(fig, 'Position', [100, 100, 500, floor(500*height/width)])
    % set(fig, 'PaperPositionMode','auto')

    
    printfig(fig, strcat(output_base_name, '_side'))

    r = 0.1 * 21.885241; 
    view(-90,90)

    frac_to_right = .35; 

    axis([-r frac_to_right*r -r r -4 2])

    set(fig, 'Position', [100, 100, 500, floor(500*(1 + frac_to_right)/2)])
    % set(fig, 'PaperPositionMode','auto')

    printfig(fig, strcat(output_base_name, '_top'))

end 
