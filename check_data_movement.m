

a = 1; 
r = 1.5;
h = 2; 
N = 4; 

filter_params.a = a; 
filter_params.r = r; 
filter_params.h = h;
filter_params.N = N;

% reference and initial surfaces are the same 
R = build_reference_surface(filter_params); 

X = R; 
alpha =  1.0; % spring constants in two directions 
beta  =  1.0;
p_0   = -1.0; 

params = pack_params(X,alpha,beta,N,p_0,R);


'version with no movement, including b.c.s'
params.X 


% 'linear order on internal points'
X_linearized = linearize_internal_points(X, params); 


'moved back to 3d ordering'
params = internal_points_to_2d(X_linearized, params);  
params.X 


















