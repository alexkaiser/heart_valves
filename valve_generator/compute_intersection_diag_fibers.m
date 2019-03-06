function coords = compute_intersection_diag_fibers(X_flat, j, k, leaflet)
% 
% Computes the intersection of the rays which have coordinates j,k
% Works for diagonal fibers only 
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
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


a = leaflet.filter.a; 


% find the indices corresponding diagonal points which lie on the ray
j_left  = leaflet.N + 2 - k; 
k_right = leaflet.N + 2 - j; 

% intersection is a solution to a linear system 
A = [a + X_flat(1,j_left,k), a - X_flat(1,j,k_right); 
         X_flat(2,j_left,k),   - X_flat(2,j,k_right)]; 

rhs = [X_flat(1,j_left,k) - X_flat(1,j,k_right); 
       X_flat(2,j_left,k) - X_flat(2,j,k_right); ]; 

soln = A \ rhs; 
s = soln(1); 

coords = [-a;0]*s + (1 - s)*X_flat(:,j_left,k); 

t = soln(2); 
coords_alt = [a;0]*t + (1 - t)*X_flat(:,j,k_right);

tol = 1e2 * eps; 
if norm(coords - coords_alt) > tol 
    error('two formulas from line disagree'); 
end 

 
