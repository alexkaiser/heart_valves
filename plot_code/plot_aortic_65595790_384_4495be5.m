

data_file = '/Users/alex/Dropbox/stanford/research_stanford/aortic_working_2019/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal/aortic_65595790_384_4495be5_circ_pt15_rad_pt54_2mm_radial_4mm_circ_circ_model_basic_updated_output_semifinal_bc_data.mat'; 

dt = 0.001665; 
nframes = 1442; 
stride = 333; 
image_size = 4*[690, 1080]; 
basename = 'aortic_65595790_384_4495be5_flow'; 

plot_aortic_movie_files(data_file, stride, nframes, basename, image_size)





