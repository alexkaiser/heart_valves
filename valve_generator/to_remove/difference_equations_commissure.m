function F = difference_equations_commissure(params, filter_params, left)
% 
% Evaluation of the global difference equations at j,k
% Uses the commissural leaflet topology
% Requires reference configuration R 
% 
% Input
%     params          Current parameters 
%     filter_params   Cone filter paramaters 
%     left            Uses the left papillary coordinate if true, right if false 
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

[X,alpha,beta,N,p_0,R,ref_frac] = unpack_params(params); 

if left 
    papillary = [0; -filter_params.a; 0]; 
else 
    papillary = [0;  filter_params.a; 0]; 
end 

F = zeros(size(X)); 

if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end 


% always 6 pressure neighbors, which may or may not be in bounds
% relative indices of pressure here 
% numbered counter clockwise 
% ignore out of bounds indices, they are not in the pressure 
% bearing part of the surface 

% here the triangulation has three types of pressure neighbors
% conveniently there are always six 
pressure_nbrs_left   = [ 0, -1; 
                         1,  0; 
                         1,  1;
                         0,  1; 
                        -1,  0;
                        -1, -1]'; 

pressure_nbrs_center = [ 0, -1; 
                         1, -1; 
                         1,  0; 
                         0,  1; 
                        -1,  0; 
                        -1, -1]'; 
                   
pressure_nbrs_right  = [ 0, -1; 
                         1, -1; 
                         1,  0; 
                         0,  1; 
                        -1,  1; 
                        -1,  0]'; 
                   
                   
for j=1:N+2
    for k=1:((N+3)/2)
        if is_internal_commissure(j,k,N)

            % pressure term first  
            pressure_term = zeros(3,1); 
            
            if j < ((N+3)/2)
                pressure_nbrs = pressure_nbrs_left; 
            elseif j == ((N+3)/2)
                pressure_nbrs = pressure_nbrs_center; 
            else
                pressure_nbrs = pressure_nbrs_right; 
            end 
                
            
            % zero indexed loop because we are computing indices with mod n 
            for n=0:5

                j_nbr      = j + pressure_nbrs(1,mod(n  ,6)+1);
                k_nbr      = k + pressure_nbrs(2,mod(n  ,6)+1);
                j_nbr_next = j + pressure_nbrs(1,mod(n+1,6)+1);
                k_nbr_next = k + pressure_nbrs(2,mod(n+1,6)+1);
                
                % fprintf('idx = %d,\t(j,k)=(%d,%d),\t(j_nbr,k_nbr)=(%d,%d),\t(j_nbr_next,k_nbr_next)=(%d,%d)\n', n, j,k, j_nbr,k_nbr, j_nbr_next,k_nbr_next); 
                
                % if any index is zero, then 
                % the pressure term does not include this value
                if j_nbr_next && k_nbr_next && j_nbr && k_nbr
                    pressure_term = pressure_term + (p_0/6) * cross(X(:,j_nbr,k_nbr) - X(:,j,k), X(:,j_nbr_next,k_nbr_next) - X(:,j,k));                     
                end 
                
            end 
            
            u_tangent_term = zeros(3,1); 
            for j_nbr = [j-1,j+1]
                    
                k_nbr = k; 

                if j_nbr == 0
                    error('no internal point should have an out of bounds j index on commisural leaflet'); 
                else 
                    X_nbr = X(:,j_nbr,k_nbr); 
                    R_nbr = R(:,j_nbr,k_nbr); 
                end 
                
                u_tangent_term = u_tangent_term + tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,alpha,ref_frac) * (X_nbr - X(:,j,k));

            end 
            
            v_tangent_term = zeros(3,1); 
            for k_nbr = [k-1,k+1]
                    
                j_nbr = j; 

                if k_nbr == 0
                    X_nbr = papillary;
                    R_nbr = papillary;
                else 
                    X_nbr = X(:,j_nbr,k_nbr); 
                    R_nbr = R(:,j_nbr,k_nbr); 
                end
                
                v_tangent_term = v_tangent_term + tension_linear(X(:,j,k),X_nbr,R(:,j,k),R_nbr,beta,ref_frac) * (X_nbr - X(:,j,k));
            end 

            F(:,j,k) = pressure_term + u_tangent_term + v_tangent_term;
            
        end
    end 
end 
    


