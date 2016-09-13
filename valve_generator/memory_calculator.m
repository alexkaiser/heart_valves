function [] = memory_calculator(N,x_mult,y_mult,z_mult)
% 
% Prints the amount of memory in Gb 
%     used for a uniform grid NS solve. 
% First estimates the size of the four fields (u,v,w,p)
% Whole thing is about a factor of ten bigger
%     since they are stored in GMRES iterations 
% 
% Input 
%     N                        Base dimension 
%     x_mult,y_mult,z_mult     Total dimension is 
%                              N * x_mult * y_mult * z_mult
% 

total_grid = 4 * N^3 * x_mult * y_mult * z_mult; 
bytes_grid = 8 * total_grid; 
GB_grid    = bytes_grid * 1e-9; 
total_mem_low  = 40 * GB_grid; 
total_mem_high  = 80 * GB_grid; 


fprintf(1, 'Memory use for grid size (%d, %d, %d):\n', N * x_mult, N * y_mult, N * z_mult); 
fprintf(1, 'Memory for one grid = %.1f Gb\n', GB_grid); 
fprintf(1, 'Total memory est (40-80 fields) = %.1f - %.1f Gb\n', total_mem_low, total_mem_high); 









