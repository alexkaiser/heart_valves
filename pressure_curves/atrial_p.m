
curve_spec = [ -8.766 -0.108 -16.498 -0.198 -22.344 -0.198 -24.104 0 -31.518 -9.888 -35.226 -9.888 -5.56 0 -66.744 15.452 -90.228 15.452 -27.368 0 -32.136 -13.596 -39.552 -13.596 -8.034 0 -12.362 14.028 -20.238 14.028 -10.044 0 -16.224 -20.826 -37.856 -20.826 -15.138 0 -25.122 11.854 -44.152 14.004 -19.03 2.15 -85.606 0.828 -110.034 0.828 -24.102 0 -31.518 -9.888 -35.226 -9.888 -5.562 0 -66.746 15.452 -90.23 15.452 -27.368 0 -32.136 -13.596 -39.552 -13.596 -8.034 0 -12.36 14.028 -20.238 14.028 -10.044 0 -16.224 -20.826 -37.856 -20.826 -18.31 0 -40.946 18.362 -62.576 12.798];

% points 5 6 are the accual coordinates 

atrial_pressure_x_RELATIVE = zeros(length(curve_spec) / 6, 1); 
atrial_pressure_y_RELATIVE = zeros(length(curve_spec) / 6, 1); 

idx = 1; 

for i=5:6:length(curve_spec)
    atrial_pressure_x_RELATIVE(idx) = curve_spec(i); 
    atrial_pressure_y_RELATIVE(idx) = curve_spec(i+1); 
    idx = idx + 1; 
end 


x_pressure = zeros(length(atrial_pressure_x_RELATIVE)+1, 1); 
y_pressure = zeros(length(atrial_pressure_y_RELATIVE)+1, 1); 


x_pressure(1) = 842.33797; 
y_pressure(1) = 419.006; 

for i=1:length(atrial_pressure_x_RELATIVE)
    x_pressure(i+1) = x_pressure(i) + atrial_pressure_x_RELATIVE(i); 
    y_pressure(i+1) = y_pressure(i) + atrial_pressure_y_RELATIVE(i); 
end 


y_pressure = y_pressure - min(y_pressure); 

range = max(y_pressure); 

% range is from -2 to 10 mmHg 
range_actual = 12; 

y_pressure = (range_actual / range) * y_pressure; 

% set correct minimum 
y_pressure = y_pressure - 2; 


x_pressure = x_pressure - min(x_pressure); 

% two beats in the current diagram 
range_actual = 1.6; 
range = max(x_pressure); 

x_pressure = (range_actual / range) * x_pressure; 

figure;
plot(x_pressure, y_pressure)
title('no smoothing')

x_spline = 0:.01:1.6; 
y_spline = spline(x_pressure, y_pressure, x_spline); 

figure; 
plot(x_spline, y_spline); 
title('atrial pressure, smoothed')


<path d="m 842.33797,432.988 
curve_spec_lv = [-32.566,-0.228 -48.238,-4.746 -52.626,-12.326 -13.488,-23.3 -7.064,-117.176 -15.472,-149.298 -4.922,-18.8 -34.346,-66.154 -70.43,-55.88 -64.85,18.462 -62.418,207.032 -69.06,209.32 -8.884,3.06 -14.372,-13.95 -28.584,-13.95 -25.34,0 -38.86,15.898 -71.528,16.738 -32.668,0.84 -115.414,4.812 -122.212,-6.93 -13.49,-23.3 -7.066,-117.176 -15.474,-149.298 -4.92,-18.8 -34.346,-66.154 -70.43,-55.88 -64.85,18.462 -62.418,207.032 -69.06,209.32 -8.886,3.06 -14.372,-13.95 -28.584,-13.95 -25.338,0 -40.17,16.738 -71.846,16.738] 


