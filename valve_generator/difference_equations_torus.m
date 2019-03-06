function [F_torus F_chordae_left F_chordae_right] = difference_equations_torus(torus)
% 
% Evaluation of the global difference equations at j,k
% Requires reference configuration R 
% 
% Input
%     torus    Current parameters 
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


X_current          = torus.X; 
N                  = torus.N; 
p_0                = torus.p_0; 
alpha              = torus.alpha; 
beta               = torus.beta;  
du                 = torus.du; 
dv                 = torus.dv; 


F_torus = zeros(size(X_current)); 

% Internal torus part 
for j=1:N
    for k=1:N
        
        % subtract one before taking mod for zero based index 
        j_minus = mod(j-1-1, N) + 1;
        j_plus  = mod(j+1-1, N) + 1;
        k_minus = mod(k-1-1, N) + 1;
        k_plus  = mod(k+1-1, N) + 1;
        
        X = X_current(:,j,k); 

        F_tmp = zeros(3,1);

        % pressure term first  
        if p_0 ~= 0
            F_tmp = F_tmp + (p_0 / (4*du*dv)) * cross(X_current(:,j_plus,k) - X_current(:,j_minus,k), X_current(:,j,k_plus) - X_current(:,j,k_minus));                     
        end 

        % u type fibers 
        for j_nbr = [j_minus,j_plus]
            
            k_nbr = k; 
            X_nbr = X_current(:,j_nbr,k_nbr); 

            F_tmp = F_tmp + alpha/du * (X_nbr-X)/norm(X_nbr-X); 

        end 

        % v type fibers 
        for k_nbr = [k_minus,k_plus]            

            j_nbr = j; 
            X_nbr = X_current(:,j_nbr,k_nbr); 

            F_tmp = F_tmp + beta/dv * (X_nbr-X)/norm(X_nbr-X); 

        end 

        F_torus(:,j,k) = F_tmp;

    end 
end

F_chordae_left  = [];  
F_chordae_right = []; 

