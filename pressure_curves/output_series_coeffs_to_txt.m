function [] = output_series_coeffs_to_txt(a_0, a_n, b_n, n, L, file_name)

if exist('file_name', 'var')
    file = fopen(file_name, 'w'); 
else 
    file = fopen('fourier_coeffs.txt', 'w'); 
end 
fprintf(file, '%d\n', n); 
fprintf(file, '%.16e\n', L); 

fprintf(file, '%.16e\n', a_0); 

for i=1:n
    fprintf(file, '%.16e    %.16e\n', a_n(i), b_n(i)); 
end 
fclose(file); 

