function run_aortic_plots(type)

if ~exist('type', 'var')
    type = '';
end 

if strcmp(type, 'high')
    fprintf("type high")
    axis_vec_pressure = [0 2.4 -5 350]; 
    axis_vec_flow     = [0 2.4 -900 900]; 
    nframes = 1442; 
    stride = 416; 
elseif strcmp(type, 'low')
    fprintf("type high")
    axis_vec_pressure = [0 2.4 -5 90]; 
    axis_vec_flow     = [0 2.4 -500 750]; 
    nframes = 1442; 
    stride = 277; 
else    
    fprintf("type default")
    nframes = 1442; 
    stride = 333; 
    axis_vec_pressure = [0 2.4 -5 160]; 
    axis_vec_flow     = [0 2.4 -700 800];
end

addpath('~/mitral_fully_discrete/scripts')
addpath('~/copies_scripts')

if isfile('bc_data.m')
    bc_data; 
elseif isfile('../bc_data.m')
    run('../bc_data.m')
else 
    error('bc_data not found'); 
end

save('bc_data.mat', 'bc_vals')

data_file = 'bc_data.mat'; 

outdir = ''; % empty for current  

% scaling was working, but for some reason an arbitrary factor of .48 = 2250/1080 larger 
image_size = 4 * [640, 1080] * (1080/2250); 
basename = 'flow_and_pressure'; 

plot_aortic_movie_files(data_file, stride, nframes, outdir, basename, image_size, axis_vec_pressure, axis_vec_flow)





