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

total_grid = 4 * N^3 * x_mult * y_mult * z_mult; 
bytes_grid = 8 * total_grid; 
GB_grid    = bytes_grid * 1e-9; 
total_mem_low  = 40 * GB_grid; 
total_mem_high  = 80 * GB_grid; 


fprintf(1, 'Memory use for grid size (%d, %d, %d):\n', N * x_mult, N * y_mult, N * z_mult); 
fprintf(1, 'Memory for one grid = %.1f Gb\n', GB_grid); 
fprintf(1, 'Total memory est (40-80 fields) = %.1f - %.1f Gb\n', total_mem_low, total_mem_high); 









