function [params] = update_leaflet_red_black(params, filter_params, tol_local, max_it_local)
%
% Full leaflet and boundary sweep. 
% Solves the nonlinear difference equations at each component. 
% Updates params in red black ordering 
% 
% Input 
%     params             Current solution parameters 
%     tol_local          Error tolerance for local newton solve 
%     max_it_local       Max iterations on each component of the Newton solve 
%     triangle_domain    Only evalutates internal points if they are below the 
%                            line y + x = 1
% 
% Output 
%     params             Current solution parameters 
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


% Loop on internal in zero based indexing 
% Start at one, end at N

% Note that internal points are 1...N
% Zero is free edge 
% j+k = N+1 is valve ring 

for j=1:params.N
    for k=1:params.N
        if mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (params.N+2)
                params = single_eqn_newton(params,filter_params,j,k,tol_local,max_it_local); 
            end 
        end 
    end 
end

for j=1:params.N
    for k=1:params.N
        if ~mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (params.N+2)
                params = single_eqn_newton(params,filter_params,j,k,tol_local,max_it_local); 
            end
            
        end 
    end 
end


    