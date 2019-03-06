function ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix)

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

file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .05; 
n_fourier_coeffs = 1000; 
plots = false; 

[a_0_ventricle a_n_ventricle b_n_ventricle Series_ventricle] = series_and_smooth(points_one_cycle_ventricle, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_ventricle(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Ventricular pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('ventricular_pressure_yellin', suffix))

% 
% fig = figure; 
% semilogy( abs(a_n_ventricle), 'k')
% hold on 
% semilogy( abs(b_n_ventricle), ':k')
% legend('ventricle cos', 'ventricle sin')
% xlabel('n')
% ylabel('|a_n|, |b_n|')
% % title('Modulus of Fourier coefficients')
% printfig(fig, 'coefficients')



[a_0_atrium a_n_atrium b_n_atrium Series_atrium] = series_and_smooth(points_one_cycle_atrium, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_atrium_series = Series_atrium(t); 
fig = figure; 
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('atrial_pressure_yellin', suffix))



fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
hold on
plot(t, vals_atrium_series, 'k'); 
title('Atrial pressure')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('both_pressure_yellin', suffix))


% t = 0:dt:(3*cycle_length); 
% vals_ventricle_series = Series_ventricle(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Pressures')
% xlabel('t')
% ylabel('p (mmHg)')
% hold on 
% vals_atrium_series = Series_atrium(t); 
% plot(t, vals_atrium_series, '--k');
% axis([0 2.4 -10 140])
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% legend('Left Ventricle', 'Left Atrium', 'Location', 'NorthWest')
% printfig(fig, 'both_pressure_yellin_three_cycles')
% 



fig = figure; 
p_diff = vals_atrium_series - vals_ventricle_series; 
plot(t, p_diff , 'k'); 
hold on 
plot(t, 0*p_diff , 'k--'); 
title('Pressure difference')
xlabel('t')
ylabel('p (mmHg)')
set(fig, 'Position', [100, 100, 1000, 500])
set(fig,'PaperPositionMode','auto')
printfig(fig, strcat('pressure_diff_yellin', suffix))



% t = 0:dt:(3*cycle_length); 
% vals_atrium_series = Series_atrium(t); 
% fig = figure; 
% plot(t, vals_atrium_series, 'k'); 
% title('Atrial pressure')
% xlabel('t')
% ylabel('p (mmHg)')
% printfig(fig, 'atrial_pressure_yellin_three_cycles')
% 
% fig = figure; 
% semilogy( abs(a_n), 'k')
% hold on 
% semilogy( abs(b_n), ':k')
% legend('atrium cos', 'atrium sin')
% xlabel('n')
% ylabel('|a_n|, |b_n|')
% title('Modulus of Fourier coefficients')



save(strcat('series_data_yellin', suffix)); 


n_coeffs_to_output = 600; 

n = n_coeffs_to_output; 

a_0 = a_0_atrium - a_0_ventricle; 
a_n = a_n_atrium - a_n_ventricle; 
b_n = b_n_atrium - b_n_ventricle; 

a_n = a_n(1:n);
b_n = b_n(1:n);
series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/cycle_length) * (1:n) .* t)' + ...  
                                 b_n .* sin((2*pi/cycle_length) * (1:n) .* t)' );   

Series_truncated = @(t) arrayfun(series_no_array, t); 

t = 0:dt:cycle_length; 
vals_ventricle_series = Series_truncated(t); 
fig = figure; 
plot(t, vals_ventricle_series, 'k'); 
title('Pressure diff, truncated series')
xlabel('t')
ylabel('p (mmHg)')

fig = figure; 
semilogy( abs(a_n), 'k')
hold on 
semilogy( abs(b_n), ':k')
legend('atrium cos', 'atrium sin')
xlabel('n')
ylabel('|a_n|, |b_n|')
title('Modulus of Fourier coefficients')


output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, file_name); 














