% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 16; 

% Show some output 
plots = false; 

% Initialize structures 
% Many parameters are in this script 

repulsive_potential = true; 


torus = initialize_torus_data_structures(N, repulsive_potential); 


% Can use a scalar pressure 
% Or a range for continuation 
p_range = torus.p_0 .* [0:.1:1]; 


% torus = solve_valve(torus, p_range); 


