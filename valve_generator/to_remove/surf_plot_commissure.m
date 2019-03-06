function fig = surf_plot_commissure(params, filter_params, left, fig)
% 
% Plots the surface and connective chordae
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

X_copy = params.X; 
N      = params.N; 

if left 
    papillary = [0; -filter_params.a; 0]; 
else 
    papillary = [0;  filter_params.a; 0]; 
end 

if mod(N,2) ~= 1
    error('must use odd N for commisural leaflet')
end 

% NaN mask in the copy 

% fill in the 3d array
% loop through N+1 to include the ring points here 
for j=1:N+2
    for k=1:((N+3)/2)
        % out of the domain? 
        % apply a NaN mask for plotting 
        if ~in_domain_with_bc(j,k,N)
            X_copy(:,j,k) = NaN; 
        end

    end 
end

% open the figure if not passed in 
if nargin < 4
    fig = figure; 
end 

x_component = squeeze(X_copy(1,:,:)); 
y_component = squeeze(X_copy(2,:,:)); 
z_component = squeeze(X_copy(3,:,:)); 

width = 1.5; 
surf(x_component, y_component, z_component, 'LineWidth',width);

axis equal 
axis auto 
hold on 

% clean up the boundaries which are ignored by surf 
string_x = [X_copy(1,1,1), X_copy(1,2,1)];
string_y = [X_copy(2,1,1), X_copy(2,2,1)];
string_z = [X_copy(3,1,1), X_copy(3,2,1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

string_x = [X_copy(1,N+1,1), X_copy(1,N+2,1)];
string_y = [X_copy(2,N+1,1), X_copy(2,N+2,1)];
string_z = [X_copy(3,N+1,1), X_copy(3,N+2,1)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 

string_x = [X_copy(1,(N+3)/2,(N+3)/2 - 1), X_copy(1,(N+3)/2,(N+3)/2)];
string_y = [X_copy(2,(N+3)/2,(N+3)/2 - 1), X_copy(2,(N+3)/2,(N+3)/2)];
string_z = [X_copy(3,(N+3)/2,(N+3)/2 - 1), X_copy(3,(N+3)/2,(N+3)/2)];
plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 


% add chordae as line segments 
k = 1; 
for j=1:(N+2)
    if is_internal_commissure(j,k,N)
        string_x = [papillary(1), X_copy(1,j,k)];
        string_y = [papillary(2), X_copy(2,j,k)];
        string_z = [papillary(3), X_copy(3,j,k)];
        plot3(string_x, string_y, string_z, 'k', 'LineWidth',width); 
    end 
end 















