
% load series_data_1e-6_cleaned_110_ventricular_max

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














