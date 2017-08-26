
% load series_data_1e-6_cleaned_110_ventricular_max

load series_data_dt_1e-6_110_ventricular_max_pressure

t = 0:.0001:(2*true_cycle_length); 
vals_atrium_series    = Series_atrium(t); 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_atrium_series, '--k'); 
hold on; 
plot(t, vals_ventricle_series, 'k'); 
legend('atrial pressure', 'ventricular pressure', 'location', 'NorthWest'); 
title('Driving pressures')
xlabel('t')
ylabel('p (mmHg)')

fig = figure; 
plot(t, vals_atrium_series - vals_ventricle_series); 
hold on 
plot(t, 0*(vals_atrium_series - vals_ventricle_series)); 

title('Pressure difference, atrium positive')
ylabel('p (mmHg)'); 
xlabel('t'); 

printfig(fig, 'difference')

figure; 

p = vals_atrium_series - vals_ventricle_series; 

p_diastole = (p >= 0) .* p;
plot(t, p_diastole)
title('pressure, positive (diastolic) only')

figure; 
frac = 1 - (p >= 0) .* (p / max(p)).^(1/2); 
plot(t, frac); 
title('fraction to systolic position'); 

min_p = min(p)
max_p = max(p)














