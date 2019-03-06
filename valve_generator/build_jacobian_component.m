function J = build_jacobian_component(params, filter_params,j,k)
% 
% Builds the Jacobian for the current index and parameter values 
% 
% Input 
%      params    Current parameter values
%      j,k       Indices, must be internal to the arrays
% 
% Output 
%      J         Jacobian of difference equations 

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

    left_papillary  = [0; -filter_params.a; 0]; 
    right_papillary = [0;  filter_params.a; 0]; 

%    F = zeros(size(X)); 

    % total internal points in triangular domain 
    % total_internal = 3*N*(N+1)/2; 
    
    J = zeros(3,3); 

    
    % always 6 pressure neighbors, which may or may not be in bounds
    % relative indices of pressure here 
    % numbered counter clockwise 
    % ignore out of bounds indices, they are not in the pressure 
    % bearing part of the surface 
    pressure_nbrs = [ 0, -1; 
                      1, -1; 
                      1,  0; 
                      0,  1; 
                     -1,  1; 
                     -1,  0]';   
                

    % in the triangle?
    if (j+k) < (N+2)

        % vertical offset does not change while differentiating this equation 
        % range_current = linear_index_offset(j,k,N) + (1:3); 


        % pressure portion 
        % zero indexed loop because we are computing indices with mod n 
        for n=0:5

            j_nbr      = j + pressure_nbrs(1,mod(n  ,6)+1); 
            k_nbr      = k + pressure_nbrs(2,mod(n  ,6)+1); 
            j_nbr_next = j + pressure_nbrs(1,mod(n+1,6)+1); 
            k_nbr_next = k + pressure_nbrs(2,mod(n+1,6)+1);

            % if any index is zero, then 
            % the pressure term does not include this triangle
            if j_nbr_next && k_nbr_next && j_nbr && k_nbr

                % Current has two terms from a product rule 
                J = J - (p_0/6) * cross_matrix(X(:,j_nbr     ,k_nbr     ) - X(:,j,k)) ...
                      + (p_0/6) * cross_matrix(X(:,j_nbr_next,k_nbr_next) - X(:,j,k));  
                      

                % 
                % No neighbors for these terms 
                % 
            end
        end 




        % u tension terms 
        for j_nbr = [j-1,j+1]

            k_nbr = k; 

            if j_nbr == 0
                X_nbr = left_papillary;
                R_nbr = left_papillary;
            else 
                X_nbr = X(:,j_nbr,k_nbr); 
                R_nbr = R(:,j_nbr,k_nbr); 
            end 

            J_tension = tension_over_norm_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,alpha,ref_frac); 

            J = J + J_tension; 

        end 


        % v tension terms 
        for k_nbr = [k-1,k+1]

            j_nbr = j; 

            if k_nbr == 0
                X_nbr = right_papillary;
                R_nbr = right_papillary;
            else 
                X_nbr = X(:,j_nbr,k_nbr); 
                R_nbr = R(:,j_nbr,k_nbr); 
            end 

            J_tension = tension_over_norm_jacobian(X(:,j,k),X_nbr,R(:,j,k),R_nbr,beta,ref_frac); 
            
            J = J + J_tension; 
 
        end 
    end
end 

