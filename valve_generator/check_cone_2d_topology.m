



a = 1; 
r = 1.5;
h = 2; 
N = 32; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;

% reference and initial surfaces are the same 
R = build_reference_surface(filter_params); 

X = R; 
alpha =  1.0; % spring constants in two directions 
beta  =  1.0;
p_0   = -1.0; % check signs on this 
ref_frac = 0.8; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac); 

fig = surf_plot(params, filter_params); 




















