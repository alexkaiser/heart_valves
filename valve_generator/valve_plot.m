function fig = valve_plot(valve, fig)
% Plots full valve data structure 

% open the figure if not passed in 
if ~exist('fig', 'var')
    fig = figure; 
end 

fig = surf_plot(valve.anterior, fig); 
hold on; 

fig = surf_plot(valve.posterior, fig); 

if isfield(valve, 'comm_left')
    fig = surf_plot(valve.comm_left, fig); 
end 

if isfield(valve, 'comm_right')
    fig = surf_plot(valve.comm_right, fig); 
end 


