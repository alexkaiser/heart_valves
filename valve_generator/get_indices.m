function [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr, k_nbr)
%
% Returns whether neighbor is a valid point, 
% If valid retuns neighbor indicies 
% and incides for spring constants 
%
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

j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 
is_bc       = leaflet.is_bc;

if isfield(leaflet, 'periodic_j')
    periodic_j = leaflet.periodic_j; 
else
    periodic_j = zeros(k_max,1); 
end

% j spring is minimum, unless zero in which case gets a periodic wrap 
j_spr = min(j, j_nbr); 
if j_spr == 0 
    j_spr = j_max; 
end 

% j_nbr may need periodic reduction 
j_nbr = get_j_nbr(j_nbr, k, periodic_j, j_max); 

% k_nbr is always an identity operation, no periodicity in this direction 
% only include for consistency 
% k_nbr = k_nbr; 

k_spr = min(k, k_nbr);

% neighbor must be valid 
if (j_nbr > 0) && (k_nbr > 0) && ...
   (j_nbr <= j_max) && (k_nbr <= k_max) && ...
   (is_internal(j_nbr,k_nbr) || is_bc(j_nbr,k_nbr))
    valid = true; 
else
    valid = false; 
end 

