function [a_0 a_n b_n Series] = fourier_series(x, y, L ,n)
%
% Fourier coefficients on [0,L]
% Simple trapezoidal quadrature 
% Allows for non-uniform spacing in x 
%
% Input: 
%    x         x inputs, must be sorted, need not be uniform 
%    y         Corresponding y values 
%    L         Interval length 
%    n         Number of coefficients to calculate 
% 
% Output: 
%    a_0       Zero cosine coefficient 
%    a_n       Cosine coefficients, does not include zero coeff 
%    b_n       Sine coefficients 
%    Series    Function handle for the series 
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


    a_n = zeros(n,1); 
    b_n = zeros(n,1); 

    c = @(t) 1/L * ones(size(t)); 
    a_0 = trap_intergral(x,y,c,L);

    for j = 1:n        
        c = @(t) (2/L) * cos((j*2*pi/L) * t); 
        a_n(j) = trap_intergral(x,y,c,L); 
        
        s = @(t) (2/L) * sin((j*2*pi/L) * t); 
        b_n(j) = trap_intergral(x,y,s,L);         
    end 
 
    series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/L) * (1:n) .* t)' + ...  
                                     b_n .* sin((2*pi/L) * (1:n) .* t)' );   

    Series = @(t) arrayfun(series_no_array, t); 

                        
%     for j = 1:n        
%         c = @(t) cos((j*2*pi/L) * t); 
%         s = @(t) sin((j*2*pi/L) * t); 
%         
%         Series = @(t) Series(t) + a_n(j) * c(t) + b_n(j) * s(t); 
%     end 
    
    
end 

function val = trap_intergral(x,y,f,L)
% simple trap rule integral 
% does not assume that x points are uniformly placed 

    N = length(x); 
    
    % mesh spacing, throws out first value which is junk
    dx = x(2:N) - x(1:N-1); 
    
    if ~all(dx > 0)
        error('Some intervals have zero or negative length'); 
    end 
    
    f_x = f(x); 
        
    val = 0.5 * sum(dx .* (y(2:N) .* f_x(2:N) + y(1:N-1) .* f_x(1:N-1))); 

    % last point periodic 
    val = val + 0.5 * ((x(N)-L) - x(1)) * (y(N)*f(x(N)-L) + y(1)*f(x(1))); 
end 














