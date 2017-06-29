function [] = output_valve_ring(valve)


X     = valve.leaflets(1).X; 
k_max = valve.leaflets.k_max; 

x_ring = squeeze(X(1,:,k_max)); 
y_ring = squeeze(X(2,:,k_max)); 

fig = figure; 
plot(x_ring, y_ring, 'k')
hold on 

% patch periodic
x_patch = [x_ring(end), x_ring(1)]; 
y_patch = [y_ring(end), y_ring(1)]; 

plot(x_patch, y_patch, 'k'); 

buffer = .3; 

axis equal; 
axis([min(x_ring)-buffer max(x_ring)+buffer min(y_ring)-buffer max(y_ring)+buffer])

xlabel('cm')
ylabel('cm')
title('Valve ring'); 

set(gcf,'color',[1 1 1])
printfig(fig, 'valve_ring_2d.eps'); 

