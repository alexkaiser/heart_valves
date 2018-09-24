
% flux_plot_IBAMR_sorted; 

times = data(:,1); 
flux = -data(:,2); 

max_abs = max(abs(flux)); 
% max_abs = 200; 
flux_scaled = flux / max_abs; 
% sound(flux_ring, sample_freq)

masking_systolic = true; 
if masking_systolic 
    % physical time of the simulation 
    t_phys = 2.4; 

    % dt is loaded from the file 
    min_t_on = 0.45; 
    full_on = 0.5;
    min_dec = 0.65; 
    full_off = 0.70; 
    beat_time = 0.8; 

    flux_masked = mask_flux(times, flux, dt, min_t_on, full_on, min_dec, full_off, beat_time); 
    figure; 
    plot(times, flux_masked); 
    title('flux masked')
    flux = flux_masked; 
end 


ode_high_pass_filtering = false; 
if ode_high_pass_filtering  
    % dt is loaded from the file 
    tau = 1e-1 / 2*pi; 
    flux_filtered = filter_sound(flux, dt, tau); 
    fig = figure; 
    plot(times, flux_filtered); 
    title('filtered flux')

    min_time_for_flux = .4; 
    min_step_for_flux = floor(min_time_for_flux/dt); 

    max_flux_filtered = max(abs(flux_filtered(min_step_for_flux:end)));

    mask = ones(size(flux_filtered)); 
    mask(1:min_step_for_flux) = linspace(0,1,min_step_for_flux); 

    flux_filtered_normalized = mask .* (flux_filtered / max_flux_filtered); 
    figure; 
    plot(times, flux_filtered_normalized); 
    title('filtered flux (normalized)')
    
    flux = flux_filtered_normalized; 
    
end 


audio_rate = 44100; 

sample_rate_good = 1/audio_rate; 
sample_points_desired = 0:sample_rate_good:t_phys; 

flux_spline = spline(times, flux, sample_points_desired); 
flux_spline_normalized = flux_spline / max(abs(flux_spline)); 

fig = figure; 
plot(sample_points_desired, flux_spline_normalized)
title('flux, as spline, to pass to sound')


sound(flux_spline_normalized, audio_rate)
audiowrite('heart_sounds_filtered.wav',flux_spline_normalized,audio_rate)




% flux_spline = spline(times, flux, sample_points_desired); 
% 
% fig = figure; 
% plot(sample_points_desired, flux_spline)
% title('flux, as spline')
% 
% flux_spline_normalized = flux_spline / max(abs(flux_spline)); 
% 
% % sound(flux_spline_normalized, audio_rate)
% audiowrite('heart_sounds.wav',flux_spline_normalized,audio_rate)
