function leaflet = get_util_arrays_bead_slip(leaflet, valve)
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

j_max                     = leaflet.j_max; 
k_min                     = leaflet.k_min; 
k_max                     = leaflet.k_max; 
ring_k_idx                = leaflet.ring_k_idx; 

N_anterior                = valve.N_anterior; 
N_posterior               = valve.N_posterior;
commissural_leaflets      = valve.commissural_leaflets; 
if commissural_leaflets 
    N_commissure          = valve.N_commissure;  
else 
    N_commissure          = 0; 
end 

% data management stuff 
is_internal               = zeros(j_max, k_max); 
is_bc                     = zeros(j_max, k_max); 
linear_idx_offset         = zeros(j_max, k_max); 
point_idx_with_bc         = zeros(j_max, k_max); 


if commissural_leaflets 
    j_range_anterior   = (1:N_anterior); 
    j_range_right_comm = (1:N_commissure) + max(j_range_anterior); 
    j_range_posterior  = (1:N_posterior)  + max(j_range_right_comm); 
    j_range_left_comm  = (1:N_commissure) + max(j_range_posterior); 
    
    indices = [j_range_anterior, j_range_right_comm, j_range_posterior, j_range_left_comm]; 
    if ~all(indices == (1:j_max))
        error('Inconsistency in indices'); 
    end 
        
else 
    j_range_anterior   = (1:N_anterior); 
    j_range_right_comm = [];
    j_range_posterior  = (1:N_posterior)  + max(j_range_anterior); 
    j_range_left_comm  = []; 
    
    indices = [j_range_anterior, j_range_posterior]; 
    if ~all(indices == (1:j_max))
        error('Inconsistency in indices'); 
    end 
    
end 


if leaflet.radial_and_circumferential 
    
    % radial and circumferential fiber layout 
    
    % valve ring at (j,ring_k_idx(j))
    for j=1:j_max 
        is_bc(j,ring_k_idx(j)) = true; 
    end 
    
    for j=1:j_max 
        
        % minimum k idx is variable 
        k = k_min(j);
        
        % keep going until hitting ring, which is always a bc 
        while ~is_bc(j,k)
            is_internal(j,k) = true; 
            k = k+1; 
        end 
    end 

else 
    error('diag fibers not implemented with bead slip')
end 


% Indices for Jacobian building 
count = 0; 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)
            linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

% Indices for spring attachments 
count = 0;
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k) || is_bc(j,k)
            point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
    end 
end


leaflet.is_internal           = is_internal;
leaflet.is_bc                 = is_bc;
leaflet.linear_idx_offset     = linear_idx_offset;
leaflet.point_idx_with_bc     = point_idx_with_bc;

leaflet.j_range_anterior   = j_range_anterior; 
leaflet.j_range_right_comm = j_range_right_comm; 
leaflet.j_range_posterior  = j_range_posterior; 
leaflet.j_range_left_comm  = j_range_left_comm; 
