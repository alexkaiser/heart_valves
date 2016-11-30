function [a_0 a_n b_n Series] = fourier_series_uniform(x, y, L, n, dt)
%
% Fourier coefficients on [0,L]
% Simple trapezoidal quadrature 
% Requires uniform spacing in x 
% Assumes that function is periodic, so just sums the values 
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

a_n = zeros(n,1); 
b_n = zeros(n,1); 

c = @(t) 1/L * ones(size(t)); 
a_0 = dt * c(x)' * y; 

for j = 1:n        
    c = @(t) (2/L) * cos((j*2*pi/L) * t); 
    a_n(j) = dt * c(x)' * y; 

    s = @(t) (2/L) * sin((j*2*pi/L) * t); 
    b_n(j) = dt * s(x)' * y; 
end 

series_no_array = @(t) a_0 + sum(a_n .* cos((2*pi/L) * (1:n) .* t)' + ...  
                                 b_n .* sin((2*pi/L) * (1:n) .* t)' );   

Series = @(t) arrayfun(series_no_array, t); 


