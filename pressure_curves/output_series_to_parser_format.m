function [] = output_series_to_parser_format(a_0_atrium, a_n_atrium, b_n_atrium, a_0_ventricle, a_n_ventricle, b_n_ventricle, n, L)

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

file = fopen('pressure_series.txt', 'w'); 


str = 'gcoef_function_4 = "-1.0*MMHG_TO_CGS*('; 
str = strcat(str, main_function_string(a_0_ventricle, a_n_ventricle, b_n_ventricle, n, L));  
str = strcat(str, ')"\n'); 

fprintf(file, str); 

fprintf(file,'\n'); 


str = 'gcoef_function_5 = "-1.0*MMHG_TO_CGS*('; 
str = strcat(str, main_function_string(a_0_atrium, a_n_atrium, b_n_atrium, n, L));  
str = strcat(str, ')"\n'); 
fprintf(file, str); 


end 




function str = main_function_string(a_0,a_n,b_n,n,L) 
    % outputs fourier series in mu parser compatible format 
    str = sprintf('%.14f ', a_0); 

    for j=1:n
        str = strcat(str, sprintf( ' + %.14f*cos(%.14f*t) + %.14f*sin(%.14f*t)', a_n(j), 2*pi*j/L, b_n(j), 2*pi*j/L)); 
    end 
end 

































