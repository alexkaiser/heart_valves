function val = cos_bump(t, angle_total, r_dip)

% reduce mod 2 pi 
t = mod(t,2*pi); 

% center on zero 
if t >= pi
    t = t - 2*pi; 
end 
    
if abs(t) < (angle_total/2) 
    val = -r_dip * cos(t * pi/angle_total)^2; 
else 
    val = 0; 
end 
    
