


% quick hack plot of atrial compliance changes 

C_LA = 1.6; 


center = .47; 
width  = .14; 
radius = width/2; 

bump = @(x) (abs(x - center) < radius) .* cos(pi * (x - center)/width).^2; 

t = 0:.001:.8; 

C_la_time_dep = @(x) C_LA * (1 - bump(x)); 


plot(t, C_la_time_dep(t))















