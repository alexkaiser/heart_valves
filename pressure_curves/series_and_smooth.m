function [a_0 a_n b_n Series times linear_interp_vals_one_cycle] = series_and_smooth(points_one_cycle, dt, bump_radius, n_fourier_coeffs, plots)
% 
% Takes data, computes piecewise linear interpolant, 
% smooths with convolution with cosine squared bump,
% computes and outputs Fourier coefficients and function handle for series 
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

debug = false; 


cycle_length = points_one_cycle(end, 1); 

times = (0:dt:(cycle_length-dt))'; 
n_times = length(times); 

% piecewise linear interpolant 
vals = interp1(points_one_cycle(:,1), points_one_cycle(:,2), times); 

linear_interp_vals_one_cycle = vals; 

if plots 
    fig = figure; 
    plot(times, vals); 
    title('Piecewise linear interpolant'); 
end 

times_three_cycle = [times; times+cycle_length; times+(2*cycle_length)]; 
vals_three_cycle = [vals; vals; vals]; 

% fig = figure; 
% plot(times_three_cycle, vals_three_cycle)
% title('three cycle')

% scaling factor that gives appropriate radius 
h = (2/pi) * bump_radius; 

% this should integrate to one 
cos_bump = @(x) (abs(x) <= bump_radius) .* (1/h) .* (2/pi) .* cos(x/h).^2; 

% want this mesh to be aligned with the previous mesh
% shift by scalar to be approx centered at zero 
mesh_bump = times_three_cycle - times_three_cycle(floor(length(times_three_cycle)/2)); 
bump_vals = cos_bump(mesh_bump); 

if debug 
    fig = figure; 
    plot(mesh_bump, bump_vals); 
    title('bump'); 

    approx_integral = dt * sum(bump_vals)  
    integral_by_quad = quad(cos_bump, -1,1) 
end 


smoothed = dt*conv(vals_three_cycle, bump_vals, 'same'); 

smoothed_one_cycle = smoothed( (n_times+1) : (2*n_times)); 

if plots
    fig = figure; 
    plot(times, smoothed_one_cycle ); 
    title('After convolution, one cycle')
end 

[a_0 a_n b_n Series] = fourier_series_uniform(times, smoothed_one_cycle, cycle_length, n_fourier_coeffs, dt); 


