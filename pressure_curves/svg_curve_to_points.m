function [x y] = svg_curve_to_points(curve_spec, x_0, y_0, min_x, max_x, min_y, max_y)
% 
% Takes coordinates from an svg curve and returns x,y pairs 
%  
% The 5th and 6th (and up in strides of 6) have the coordinates 
% 
% Initial coordinates must be provided 
% Movement is relative 
% 




x_RELATIVE = zeros(length(curve_spec) / 6, 1); 
y_RELATIVE = zeros(length(curve_spec) / 6, 1); 

if mod(length(curve_spec,6)) ~= 0
    error('must give 6 points per coordinate')
end

% points 5 6 are the accual coordinates 
idx = 1;
for i=5:6:length(curve_spec)
    x_RELATIVE(idx) = curve_spec(i); 
    y_RELATIVE(idx) = curve_spec(i+1); 
    idx = idx + 1; 
end 

x = zeros(length(x_RELATIVE)+1, 1); 
y = zeros(length(y_RELATIVE)+1, 1); 


x(1) = x_0; 
y(1) = y_0; 

for i=1:length(x_RELATIVE)
    x(i+1) = x(i) + x_RELATIVE(i); 
    y(i+1) = y(i) + y_RELATIVE(i); 
end 

% zero the minimum, scale the array, fix the real minimum
x = x - min_x; 
range = max(x); 
range_actual = max(x) - min(x); 
x = (range_actual / range) * x;
x = x + min_x; 

y = y - min_y; 
range = max(y); 
range_actual = max(y) - min(y); 
y = (range_actual / range) * y; 
y = y + min_y; 











