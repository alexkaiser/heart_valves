function flux_masked = mask_flux(times, flux, dt, min_t_on, full_on, min_dec, full_off, beat_time)
% 
% masks flux with piecewise linear function 
% function is zero, then linearly to one, then linearly back to zero 
% mask is repeated periodically 
%
% times - times of flux 
% flux - values 
% min_t_on - linear increase from zero to one here  
% full_on - mask is one starting at this time 
% min_dec - mask goes linearly to zero starting here 
% full_off - mask is zero until end of cycle here 
% beat_time - time with which to repeat this 
% 

t_interp = [0, min_t_on, full_on, min_dec, full_off]; 
values   = [0,        0,       1,       1,        0]; 
if ~issorted(times)
    error('must provide a sorted list of times'); 
end 

% if we get an extra beat here, no problem 
beats = ceil(length(times) * dt / beat_time); 

for i=1:(beats-1)
    t_interp_temp = [0, min_t_on, full_on, min_dec, full_off] + (beat_time*i); 
    values_temp   = [0,        0,       1,       1,        0]; 
    t_interp = [t_interp, t_interp_temp]; 
    values   = [values, values_temp]; 
end

if t_interp(end) ~= max(times)
    t_interp = [t_interp, max(times)]; 
    values   = [values,           0]; 
end 

mask = interp1(t_interp, values, times); 

flux_masked = flux .* mask; 

