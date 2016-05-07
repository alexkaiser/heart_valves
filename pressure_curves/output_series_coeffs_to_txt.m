function [] = output_series_coeffs_to_txt(a_0_atrium, a_n_atrium, b_n_atrium, a_0_ventricle, a_n_ventricle, b_n_ventricle, n, L)

file = fopen('fourier_coeffs.txt', 'w'); 

fprintf(file, '%d\n', n); 
fprintf(file, '%.16e\n', L); 

a_0 = a_0_atrium - a_0_ventricle; 
a_n = a_n_atrium - a_n_ventricle; 
b_n = b_n_atrium - b_n_ventricle; 

fprintf(file, '%.16e\n', a_0); 

for i=1:n
    fprintf(file, '%.16e    %.16e\n', a_n(i), b_n(i)); 
end 
fclose(file); 

