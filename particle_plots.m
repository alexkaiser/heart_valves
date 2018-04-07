

data_dir = '/Users/alex/mitral_fully_discrete/plots_and_session_files/mitral_1647284_512_git_5855699_two_leaflet_mean_reverting'

file_name = 'mitral_tree_512_001005.m'

max_velocity = 200; 

clear_fucntion_cache = false; 

if clear_fucntion_cache 
    file_name_copy = file_name; 
    clear file_name
    file_name = file_name_copy;
end 

file_name = [data_dir '/' file_name]

fig = plot_particles(file_name, max_velocity); 
printfig(fig, 'particle_test')


