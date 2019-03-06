function coords = compute_intersection_commissure(X_flat, j, k, filter_params, left)
% 
% Computes the intersection of the rays which have coordinates j,k
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
%     left            Uses left papillary if true
% 
% Output: 
%     coords          2d vector for intersection 
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


if left 
    papillary = [-filter_params.a; 0]; 
else 
    papillary = [ filter_params.a; 0]; 
end 

N = filter_params.N; 
if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end

% k_vert is the k on the b.c. above the current j 
if j <= ((N+3)/2)
    k_vert = j; 
else 
    k_vert = N + 3 - j; 
end 

% j_horiz is the j on the right boundary 
j_horiz = N + 3 - k; 

% the left boundary is always k,k so no extra term needed 

% intersection is a solution to a linear system 
A = [X_flat(1,k,k) - X_flat(1,j_horiz,k), -papillary(1) + X_flat(1,j,k_vert);  
     X_flat(2,k,k) - X_flat(2,j_horiz,k), -papillary(2) + X_flat(2,j,k_vert)];   

rhs = [-X_flat(1,j_horiz,k) + X_flat(1,j,k_vert); 
       -X_flat(2,j_horiz,k) + X_flat(2,j,k_vert)]; 

soln = A \ rhs; 
s = soln(1); 

coords = s*X_flat(:,k,k) + (1 - s)*X_flat(:,j_horiz,k); 

t = soln(2); 
coords_alt = t*papillary + (1 - t)*X_flat(:,j,k_vert);

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 



