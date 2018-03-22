function cmap = make_colormap(len)
% 
% returns a custom color map of length len 
% 
% 
% 
% 

hue = linspace(4/6,0,len)'; 

min_full_saturation = floor(len * 3/4); 

saturation = linspace(0,1,min_full_saturation)'; 
saturation = [saturation; ones(len-min_full_saturation,1);]; 

max_full_value = floor(len * 3/4); 
value = ones(max_full_value, 1); 
value = [value; linspace(1,0,len-max_full_value)']; 

hsv = [hue, saturation, value]; 

cmap = hsv2rgb(hsv); 
