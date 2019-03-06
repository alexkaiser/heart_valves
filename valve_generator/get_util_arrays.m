function [is_internal is_bc linear_idx_offset point_idx_with_bc] = get_util_arrays(leaflet)
% 
% Returns three arrays with information about the geometry 
% 
% Output: 
%     is_internal          Boolean, true if 
%     is_bc                Point is a boundary condition that is fixed 
%     linear_idx_offset    In Jacobian, linear_idx_offset(j,k) + 1:3
%                          are the indices for the vector X(:,j,k)
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

N                       = leaflet.N; 
j_max                   = leaflet.j_max; 
k_max                   = leaflet.k_max; 

if isfield(leaflet, 'trapezoidal_flat_points')
    trapezoidal_flat_points = leaflet.trapezoidal_flat_points; 
else 
    trapezoidal_flat_points = 0; 
end 

is_internal       = zeros(j_max, k_max); 
is_bc             = zeros(j_max, k_max); 
linear_idx_offset = zeros(j_max, k_max); 
point_idx_with_bc = zeros(j_max, k_max); 


if leaflet.radial_and_circumferential 
    
    % radial and circumferential fiber layout 
    
    % valve ring at k_max
    k=k_max; 
    for j=1:j_max 
        is_bc(j,k) = true; 
    end 
    
    % loop from free edge then up in k 
    j = k_max;  
    for k=1:(k_max - 1)
        for k_tmp=k:(k_max-1)
            is_internal(j,k_tmp) = true; 
        end 
        j = j - 1; 
    end 

    for j = (k_max+1):(k_max + trapezoidal_flat_points)
        for k=1:(k_max-1)
            is_internal(j,k) = true; 
        end 
    end 
    
    j = k_max + 1 + trapezoidal_flat_points; 
    for k=1:(k_max - 1)
        for k_tmp=k:(k_max-1)
            is_internal(j,k_tmp) = true; 
        end 
        j = j + 1; 
    end 
        
else 

    % diagonal fiber layout 
        
    is_internal       = zeros(j_max, k_max); 
    is_bc             = zeros(j_max, k_max); 
    linear_idx_offset = zeros(j_max, k_max); 
    point_idx_with_bc = zeros(j_max, k_max); 
    
    k = N+1; 
    for j=1:j_max
        is_bc(j,k) = true;
        k = k-1; 
    end 
    
    for j=1:j_max
        for k=1:k_max
            % in the triangle? 
            if ((j+k) < (N+2))
                is_internal(j,k) = true;  
            end 
        end 
    end 

end 


count = 0; 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

count = 0;
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k) || is_bc(j,k)
            point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
    end 
end
