function [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k)
% 
% Returns neighbors for computing pressure  
% The point j,k should be internal (not a boundary condition)
% 
% If the point is not on the free edge, then points up and down
% for a two mesh widths wide, centered finite difference are supplied 
% 
% If the point is on the free edge, then a 
% 
% The value "m" is the coefficient of the finite difference pressure 
% This is four internal to the leaflet
% but is 1/2 if one direction uses a one mesh width stencil  
% and 1 if both directions use one mesh width stencils 
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

num_j_modified = 0; 
[valid j_plus__1] = get_indices(leaflet, j, k, j+1, k); 
if ~valid
    j_plus__1 = j; 
    num_j_modified = num_j_modified + 1; 
end 

[valid j_minus_1] = get_indices(leaflet, j, k, j-1, k); 

if ~valid
    j_minus_1 = j; 
    num_j_modified = num_j_modified + 1; 
end 

num_k_modified = 0; 
[valid tmp k_plus__1] = get_indices(leaflet, j, k, j, k+1); 
if ~valid
    k_plus__1 = k; 
    num_k_modified = num_k_modified + 1; 
end 

[valid tmp k_minus_1] = get_indices(leaflet, j, k, j, k-1); 

if ~valid
    k_minus_1 = k; 
    num_k_modified = num_k_modified + 1; 
end 


if     (num_j_modified == 0) && (num_k_modified == 0)
    m = 0.25; 
elseif (num_j_modified == 1) && (num_k_modified == 0)
    m = 0.5; 
elseif (num_j_modified == 0) && (num_k_modified == 1)
    m = 0.5; 
elseif (num_j_modified == 1) && (num_k_modified == 1)
    m = 1.0; 
else 
    error('All pressure indices must be modified zero or one'); 
end 




