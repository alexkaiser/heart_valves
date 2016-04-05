function [Series] = fourier_series_function_handle(a_0, a_n, b_n, n, L)



series_no_array = @(t) a_0 + sum(a_n(1:n) .* cos((2*pi/L) * (1:n) .* t)' + ...  
                                 b_n(1:n) .* sin((2*pi/L) * (1:n) .* t)' );   

Series = @(t) arrayfun(series_no_array, t);