
% 0 

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

t_diastole_full = .1;  
t_systole_start = .46; % .44;  
t_systole_full  = .5;  

t_cycle_length  = .8; 

times = 0:.001:t_cycle_length; 
vals = zeros(size(times)); 
derivs = zeros(size(times)); 

derivs_linear = zeros(size(times)); 


for i=1:length(times)

    current_time = times(i); 
    
    t_reduced = current_time - t_cycle_length * floor(current_time/(t_cycle_length)); 
    
    
    if (t_reduced < t_diastole_full)
        vals(i)          = (1/2) * (1 - cos(pi * t_reduced/t_diastole_full)); 
        derivs(i)        = pi/(2*t_diastole_full) * sin(pi * t_reduced/t_diastole_full); 
        derivs_linear(i) = 1/t_diastole_full; 
        
        
    elseif  (t_reduced < t_systole_start)
        
        vals(i)          = 1; 
        derivs(i)        = 0; 
        derivs_linear(i) = 0; 
   
    elseif  (t_reduced < t_systole_full)
        
        vals(i)          =  (1/2) * (cos(pi * (t_reduced - t_systole_start)/(t_systole_full - t_systole_start)) + 1); 
        derivs(i)        = -(1/2) *  sin(pi * (t_reduced - t_systole_start)/(t_systole_full - t_systole_start)) * (pi/(t_systole_full - t_systole_start)) ; 
        
        derivs_linear(i) = -1 / (t_systole_full - t_systole_start); 
        
    else
    
        vals(i)          = 0; 
        derivs(i)        = 0;  
        derivs_linear(i) = 0; 
        
    end 

end 


figure; 
plot(times, vals,'k' )
hold on 
plot(times, derivs, 'k--')

plot(times, derivs_linear, 'k-.')

legend('values', 'derivatives', 'piecewise linear deriv')



