function fig = valve_plot(valve, fig)
% Plots full valve data structure 

% open the figure if not passed in 
if ~exist('fig', 'var')
    fig = figure; 
end 


for i=1:length(valve.leaflets)
    
    fig = surf_plot(valve.leaflets(i), fig); 
    axis equal
    hold on 
    
end 

hold off 
