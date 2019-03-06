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

% Copyright (c) 2019, Alexander D. Kaiser
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

