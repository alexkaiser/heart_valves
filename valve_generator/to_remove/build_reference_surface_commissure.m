function X = build_reference_surface_commissure(filter_params, left)
%
% Builds reference coffee cone surface 
% 
% Mesh has a triangular layout 
% j==1 and k==1 are assumed to be the free edge 
% j+k == N+2 is the valve ring 
% 
% 
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

N = filter_params.N; 
r = filter_params.r; 
h = filter_params.h; 

if (filter_params.min_angle < -pi) || (filter_params.max_angle > pi) 
    error('outside allowable range of angles for current parameters'); 
end 

if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end 

X      = zeros(3,N+2,(N+3)/2); 
X_flat = zeros(2,N+2,(N+3)/2); 

mesh = linspace(filter_params.min_angle,filter_params.max_angle,N+2); 
ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

% set the valve ring 
for j=1:((N+3)/2)
    k = j; 
    X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
    X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
end 

for j=((N+3)/2 + 1):(N+2)
    k = N + 3 - j; 
    X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
    X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
end 

% fill in the 3d array 
for j=1:N+2
    for k=1:((N+3)/2)
        if is_internal_commissure(j,k,N)
            X_flat(:,j,k) = compute_intersection_commissure(X_flat, j, k, filter_params, left); 
            X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
        end

    end 
end





















