function [F_leaflet F_chordae_left F_chordae_right] = difference_equations_bead_slip_leaflet_only(leaflet)
% 
% Evaluation of the global difference equations at j,k
% Requires reference configuration R 
% 
% Input
%     leaflet    Current parameters 
%
% Output
%     F        Values of all difference equation, 3 by triangular array 
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

X_current          = leaflet.X; 
p_0                = leaflet.p_0; 
alpha              = leaflet.alpha; 
beta               = leaflet.beta;  
chordae_idx_left   = leaflet.chordae_idx_left; 
chordae_idx_right  = leaflet.chordae_idx_right;
j_max              = leaflet.j_max; 
k_max              = leaflet.k_max; 
du                 = leaflet.du; 
dv                 = leaflet.dv; 
is_internal        = leaflet.is_internal; 


F_leaflet = zeros(size(X_current)); 

S_left  = alpha * ones(k_max-1,1); 
S_right = alpha * ones(k_max-1,1);     
T       = beta  * ones(j_max,1); 

% Tensions are equalized from left to right 
S = (S_left + S_right)/2.0;

% Convert from units of force to force densities 
% S = S/dv;
% T = T/du;

% Internal leaflet part 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) && (~chordae_idx_left(j,k)) && (~chordae_idx_right(j,k))

            X = X_current(:,j,k); 

            F_tmp = zeros(3,1);

            % pressure term first  
            if p_0 ~= 0
                F_tmp = F_tmp + (p_0 / (4*du*dv)) * cross(X_current(:,j+1,k) - X_current(:,j-1,k), X_current(:,j,k+1) - X_current(:,j,k-1));                     
            end 

            % u type fibers 
            for j_nbr = [j-1,j+1]

                k_nbr = k; 
                X_nbr = X_current(:,j_nbr,k_nbr); 

                F_tmp = F_tmp + S(k)/du * (X_nbr-X)/norm(X_nbr-X); 

            end 

            % v type fibers 
            for k_nbr = [k-1,k+1]

                j_nbr = j; 
                X_nbr = X_current(:,j_nbr,k_nbr); 

                F_tmp = F_tmp + T(j)/dv * (X_nbr-X)/norm(X_nbr-X); 

            end 

            F_leaflet(:,j,k) = F_tmp;

        end
    end 
end

 

F_chordae_left  = [];  
F_chordae_right = []; 

