function [] = check_pressure_curve()

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

    MMHG_TO_CGS = 1333.22368; 
    
    % pressures, higher pressure in atrium is positive 
    diastolic_p =   10.0 * MMHG_TO_CGS; 
    systolic__p = -100.0 * MMHG_TO_CGS; 

    % times for opening, open, closing and closed  
    t_opening  = 0.08; 
    t_diastole = 0.5; 
    t_closing  = 0.02;
    t_systole  = 0.2; 

    beat_time = t_diastole + t_systole + t_opening + t_closing; 

    if beat_time ~= 0.8
        warning('total beat time is not .8')
    end 

    % slope on opening part of curve 
    a_opening = (diastolic_p - systolic__p)/t_opening; 
    a_closing = (systolic__p - diastolic_p)/t_closing; 
    
    % time at which valve is completely open 
    t_1 = diastolic_p / a_opening; 
    
    % time at which closing begins 
    t_2 = t_1 + t_diastole; 
    
    % linear offset for closing part of function 
    b = diastolic_p - a_closing * t_2; 
    
    % time at which constant systolic pressure is achieved 
    t_3 = t_2 + t_closing; 
    
    % time at which opening begins 
    t_4 = t_3 + t_systole; 
    
    % linear offset for opening part of function at end 
    c = systolic__p - a_opening * t_4; 
    
    dt = 1e-4; 
    
    times = 0:dt:2*beat_time; 
    pressure = zeros(size(times)); 
    
    for j=1:length(times)
        
        if j == length(times)
            'cool'
        end 
        
        pressure(j) = pressure_curve(times(j)); 
    end 
    
    figure; 
    plot(times, pressure); 
    title('pressure with atrial side positive, cgs units'); 

    
    
    times = t_2:dt:t_2+beat_time; 
    pressure = zeros(size(times)); 
    
    for j=1:length(times)
        
        if j == length(times)
            'cool'
        end 
        
        pressure(j) = pressure_curve(times(j)); 
    end 
    
    figure; 
    plot(times, -pressure/MMHG_TO_CGS); 
    title('pressure with ventricle side positive, start at MV close, mmHg'); 
    
    
    
    function p = pressure_curve(t)
        % piecewise linear function for pressure 
        
        t_reduced  = t - beat_time*round(t/beat_time - 0.5); 
        
        if t_reduced < t_1 
            p = a_opening * t_reduced; 
        elseif t_reduced < t_2
            p = diastolic_p; 
        elseif t_reduced < t_3
            p = a_closing * t_reduced + b; 
        elseif t_reduced < t_4
            p = systolic__p; 
        else 
            p = a_opening * t_reduced + c; 
        end 
        
    end 
       
end 






