function [] = output_series_to_parser_format(a_0_atrium, a_n_atrium, b_n_atrium, a_0_ventricle, a_n_ventricle, b_n_ventricle, n, L)

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

































