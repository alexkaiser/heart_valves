% Script to build valve 


% Size parameter
% Number of points on free edge of each leaflet 
% 
N = 32; 

% Show some output 
plots = true; 

% Initialize structures 
% Many parameters are in this script 
valve = initialize_valve_data_structures(N); 

if plots 
    fig = surf_plot(valve.posterior); 
    title('Reference configuration of posterior surface'); 
    fig = surf_plot(valve.anterior); 
    title('Reference configuration of anterior surface'); 
end
    
% Can use a scalar pressure 
% Or a range for continuation 
p_range = valve.p_0; 























