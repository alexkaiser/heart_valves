

a = 1; 
r = 1.5;
h = 2; 
N = 16 + 1; 
extra = pi/4; 
center = -pi/2; 
min_angle =  center - extra; 
max_angle =  center + extra; 
left = true; 


filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;
filter_params.min_angle = min_angle;
filter_params.max_angle = max_angle;



% reference and initial surfaces are the same 
R = build_reference_surface_commissure(filter_params, left );
    
X = R; 
alpha =  1.0; % spring constants in two directions 
beta  =  1.0;
p_0   = -1.0; % check signs on this 
ref_frac = 0.8; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

fig = surf_plot_commissure(params, filter_params,left); 
title('left commisural leaflet in reference configuration')
