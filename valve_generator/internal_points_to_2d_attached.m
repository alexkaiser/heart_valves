function valve = internal_points_to_2d_attached(v_linearized, valve)
%
%  Takes the internal values in X which are arranged in linear order
%  And places them back in the 3d vector array in leaflet
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

idx = 1; 

j_max       = valve.anterior.j_max; 
k_max       = valve.anterior.k_max; 
is_internal = valve.anterior.is_internal; 

total_internal = 3*sum(is_internal(:)); 

% k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            valve.anterior.X(:,j,k) = v_linearized(idx + (0:2)); 
            idx = idx + 3; 
        end 
    end 
end 


j_max       = valve.posterior.j_max; 
k_max       = valve.posterior.k_max; 
is_internal = valve.posterior.is_internal; 

total_internal = total_internal + 3*sum(is_internal(:)); 

% k is required to be the outer loop 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            valve.posterior.X(:,j,k) = v_linearized(idx + (0:2)); 
            idx = idx + 3; 
        end 
    end 
end 

% copy chordae if length allows 
C_left   = valve.anterior.chordae.C_left; 
C_right  = valve.anterior.chordae.C_right; 

[m N_chordae] = size(C_left); 

idx = total_internal + 1; 
for i=1:N_chordae
    C_left(:,i)  = v_linearized(idx + (0:2));  
    idx = idx + 3; 
end 

idx = total_internal + 3*N_chordae + 1; 
for i=1:N_chordae
    C_right(:,i) = v_linearized(idx + (0:2));  
    idx = idx + 3; 
end 

valve.anterior.chordae.C_left  = C_left; 
valve.anterior.chordae.C_right = C_right; 


