

beat_time = 0.8; 

t = 0:.001:(2*beat_time); 

% systole_middle = .6; 

% f = @(t) sin( (pi/beat_time) * t  - systole_middle).^100; 

half_width = .15; 

bottom_support = .48;
systole_middle = bottom_support + half_width; 


bump = @(t) 0.5*(cos(pi*t) + 1) .* (abs(t) < 1.0); 

t_reduced = @(t) t - beat_time*floor(t/beat_time); 

figure; 
plot(t,beat_time*floor(t/beat_time))

figure 
plot(t, t_reduced(t))

f = @(t) bump( ( t_reduced(t) - systole_middle) / half_width);    


figure
plot(t, f(t))









