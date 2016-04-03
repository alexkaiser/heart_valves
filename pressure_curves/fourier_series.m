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

    a_n = zeros(n,1); 
    b_n = zeros(n,1); 

    c = @(t) 1/L; 
    a_0 = trap_intergral(x,y,c);

    for j = 1:n        
        c = @(t) (2/L) * cos((j*2*pi/L) * t); 
        a_n(j) = trap_intergral(x,y,c); 
        
        s = @(t) (2/L) * sin((j*2*pi/L) * t); 
        b_n(j) = trap_intergral(x,y,s);         
    end 
 
    Series = @(t) a_0; 
    
    for j = 1:n        
        c = @(t) cos((j*2*pi/L) * t); 
        s = @(t) sin((j*2*pi/L) * t); 
        
        Series = @(t) Series(t) + a_n(j) * c(t) + b_n(j) * s(t); 
    end 
    
    
end 

function val = trap_intergral(x,y,f)
% simple trap rule integral 
% does not assume that x points are uniformly placed 

    val = 0.0; 
    N = length(x); 
    for i=1:N-1
        val = val + 0.5 * (x(i+1) - x(i)) * (y(i+1)*f(x(i+1)) + y(i)*f(x(i))); 
    end 

    % val = val + 0.5 * ((x(N)-L) - x(1)) * (y(N)*f(x(N)-L) + y(1)*f(x(1))); 
end 














