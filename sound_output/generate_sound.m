
flux_plot_IBAMR_sorted_one_plot; 

times = data(:,1); 
flux = -data(:,2); 


% stride = 278; 
% 
% t_per_output = dt * stride; 
% sample_freq = 1/t_per_output; 
% soundsc(flux_ring, sample_freq)


max_abs = max(abs(flux)); 
% max_abs = 200; 
flux_scaled = flux / max_abs; 
% sound(flux_ring, sample_freq)

% physical time of the simulation 
t_phys = 2.4; 


audio_rate = 44100; 

sample_rate_good = 1/audio_rate; 

sample_points_desired = 0:sample_rate_good:t_phys; 

flux_spline = spline(times, flux, sample_points_desired); 

fig = figure; 
plot(sample_points_desired, flux_spline)
title('flux, as spline')

flux_spline_normalized = flux_spline / max(abs(flux_spline)); 

sound(flux_spline_normalized, audio_rate)


log_filter = @(x) sign(x) .* (((abs(x) < 1) .* abs(x)) + ((abs(x) > 1) .* log(abs(x)))); 

flux_log_filtered   = log_filter(flux_spline); 
fig = figure; 
plot(sample_points_desired, flux_log_filtered)
title('flux with log filter, non-normalized');
flux_log_filtered_normalized = flux_log_filtered / max(abs(flux_log_filtered)); 
sound(flux_log_filtered_normalized, audio_rate)
fig = figure; 
plot(sample_points_desired, flux_log_filtered_normalized)
title('flux with log filter'); 

flux_log_log_filtered   = log_filter(flux_log_filtered); 
fig = figure; 
plot(sample_points_desired, flux_log_log_filtered)
title('flux with log log filter, non-normalized');
flux_log_log_filtered_normalized = flux_log_log_filtered / max(abs(flux_log_log_filtered)); 
sound(flux_log_log_filtered_normalized, audio_rate)
fig = figure; 
plot(sample_points_desired, flux_log_filtered_normalized)
title('flux with log log filter');


n = 100; 
nth_root_filter = @(x) sign(x) .* abs(x).^(1/n); 

flux_nth_root_filtered   = nth_root_filter(flux_spline_normalized); 

flux_nth_root_filtered_normalized = flux_nth_root_filtered / max(abs(flux_nth_root_filtered)); 
sound(flux_nth_root_filtered, audio_rate)
fig = figure; 
plot(sample_points_desired, flux_nth_root_filtered)
title('flux with nth root filter');


% audiowrite('heart_sounds.wav',flux_spline_normalized,audio_rate)

f = fopen('time_series.txt', 'w'); 

for i=1:length(flux_spline)
    fprintf(f, '%.14f ', flux_spline_normalized(i)); 
end 
fprintf(f, '\n'); 

fclose(f); 




