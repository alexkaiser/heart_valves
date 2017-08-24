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


% add two little ticks 
tan = [x_ring(end) - x_ring(2), y_ring(end) - y_ring(2)]; 

normal = [tan(2) -tan(1)]; 
normal = normal / norm(normal); 

len = .15; 

x_tick = [x_ring(1) + len*normal(1), x_ring(1) - len*normal(1)]; 
y_tick = [y_ring(1) + len*normal(2), y_ring(1) - len*normal(2)]; 
plot(x_tick, y_tick, 'k') 


idx_right = valve.leaflets(1).j_range_anterior(end); 
tan = [x_ring(idx_right+1) - x_ring(idx_right-1), y_ring(idx_right+1) - y_ring(idx_right-1)]; 

normal = [tan(2) -tan(1)]; 
normal = normal / norm(normal); 

x_tick = [x_ring(idx_right) + len*normal(1), x_ring(idx_right) - len*normal(1)]; 
y_tick = [y_ring(idx_right) + len*normal(2), y_ring(idx_right) - len*normal(2)]; 
plot(x_tick, y_tick, 'k') 





buffer = .3; 

axis equal; 
axis([min(x_ring)-buffer max(x_ring)+buffer min(y_ring)-buffer max(y_ring)+buffer])

xlabel('cm')
ylabel('cm')
title('Valve ring'); 

set(gcf,'color',[1 1 1])
printfig(fig, 'valve_ring_2d.eps'); 

