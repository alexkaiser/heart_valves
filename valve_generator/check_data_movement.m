

a = 1; 
r = 1.5;
h = 2; 
N = 4; 
min_angle = -pi/2; 
max_angle =  pi/2; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;
filter_params.min_angle = min_angle;
filter_params.max_angle = max_angle;

% reference and initial surfaces are the same 
R = build_reference_surface(filter_params); 

X = R; 
alpha =  1.0; % spring constants in two directions 
beta  =  1.0;
p_0   = -1.0; 
ref_frac = 0.8; 

params = pack_params(X,alpha,beta,N,p_0,R,ref_frac);


'version with no movement, including b.c.s'
params.X 


% 'linear order on internal points'
X_linearized = linearize_internal_points(X, params); 


'moved back to 3d ordering'
params = internal_points_to_2d(X_linearized, params);  
params.X 


chordae_tree = true; 

if chordae_tree
    k_0 = 1; 
    k_multiplier = 2; 
    tree_frac = 0.5; 
    params = add_chordae(params, filter_params, k_0, k_multiplier, tree_frac); 
    chordae = params.chordae;
    [m N_chordae] = size(params.chordae.C_left); 
end 

'with chordae before linearization'
X 
params.chordae.C_left
params.chordae.C_right 

'linearized with chordae'
X_and_chordae_linearized = linearize_internal_points(X, params, params.chordae.C_left, params.chordae.C_right) 


'after return to normal data structure'
params = internal_points_to_2d(X_and_chordae_linearized, params);  
params.X 
params.chordae.C_left
params.chordae.C_right 











