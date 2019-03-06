% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 

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


cycle_length = 0.8; 

% quadrature spacing 
dt = 1e-5; 

plots = true; 

base_name = 'fourier_coeffs';


basic = false; 
if basic 

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120; 
    0.75, 130; 
    cycle_length, 8]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 20; 
    0.53, 5; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = ''; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 

low = true; 
if low 

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120*.5; 
    0.75, 130*.5; 
    cycle_length, 8]; 

    depressed_low_pressure = false; 
    if depressed_low_pressure
        points_one_cycle_atrium = [0.0, 12.5; 
        0.06, 4; 
        0.40, 7; 
        0.47, 20; 
        0.53, 5; 
        0.58, 7; 
        0.7,  8; 
        cycle_length, 12.5];
    else 
        points_one_cycle_atrium = [0.0, 24.555; 
        0.06, 4; 
        0.40, 7; 
        0.47, 20; 
        0.53, 5; 
        0.58, 7; 
        0.7,  10; 
        cycle_length, 24.555]; 
    end 
    
    suffix = '_low'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 

end 


high = false; 
if high 

    points_one_cycle_ventricle = [0.0,   -4; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120 * 2; 
    0.75, 130 * 2; 
    cycle_length, -4]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 20; 
    0.53, 5; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = '_high'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 


no_atrial_systole = false; 
if no_atrial_systole

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 6; 
    0.58, 120; 
    0.75, 130; 
    cycle_length, 8]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 7; 
    0.53, 7; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = 'no_atrial_systole'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 
