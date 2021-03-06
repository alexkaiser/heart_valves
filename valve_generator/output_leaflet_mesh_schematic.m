function [] = output_leaflet_mesh_schematic(valve)

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

leaflet = valve.leaflets(1); 
N                      = leaflet.N; 
X_current              = leaflet.X; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max;
j_max_anterior         = max(leaflet.j_range_anterior); 
du                     = leaflet.du; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
n_rings_periodic       = leaflet.n_rings_periodic; 

X_schematic = nan * zeros(size(X_current)); 

plot_options = 'k-'; 
marker_size = 3; 


fig = figure; 
hold on; 


N_anterior  = j_max_anterior; 
N_posterior = N - N_anterior; 


for j=1:j_max
    for k=1:k_max
        if is_internal(j, k) || is_bc(j,k)
            
            % use the real locations in u,v from order check here 
            if (j <= j_max_anterior)
                u = .5 * (j-1) / (N_anterior-1);
            else 
                j_reduced = j - j_max_anterior;
                u = .5 + .5 * j_reduced / (N_posterior + 1);  
            end 
            
            v = 1 - (k_max - k) * du; 
            
            X_schematic(:,j,k) = [u,v,nan]; %du * [j;k;nan]; 
        end 
    end 
end 

leaflet_schematic   = leaflet; 
leaflet_schematic.X = X_schematic; 


for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) 

            % plot all connections to current location 
            % redundantly becuase who cares 

            % u type fibers 
            for j_nbr_tmp = [j-1,j+1]

                k_nbr_tmp = k; 

                [valid j_nbr k_nbr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid 
                    x_tmp(1) = X_schematic(1,j,k); 
                    x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                    y_tmp(1) = X_schematic(2,j,k); 
                    y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 

                    % lazy hack to not plot weird periodic links 
                    if abs(x_tmp(1) - x_tmp(2)) < 2*du
                        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                    end 

                end 
            end 

            % v type fibers 
            for k_nbr_tmp = [k-1,k+1]

                j_nbr_tmp = j; 

                [valid j_nbr k_nbr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

                if valid
                    x_tmp(1) = X_schematic(1,j,k); 
                    x_tmp(2) = X_schematic(1,j_nbr,k_nbr); 
                    y_tmp(1) = X_schematic(2,j,k); 
                    y_tmp(2) = X_schematic(2,j_nbr,k_nbr); 
                    plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 
                end 
            end 
        end       
    end 
end 

% more hacks, clean up periodic boundary 
j=j_max; 
for k=1:k_max
    if is_internal(j,k) 
        x_tmp(1) = X_schematic(1,j,k); 
        x_tmp(2) = X_schematic(1,j,k) + du/2; 
        y_tmp(1) = X_schematic(2,j,k); 
        y_tmp(2) = X_schematic(2,j,k); 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

        % make a little zig zag to show periodicity 
        x_tmp(1) = X_schematic(1,j,k) + du/2; 
        x_tmp(2) = X_schematic(1,j,k) + 3*du/4; 
        y_tmp(1) = X_schematic(2,j,k); 
        y_tmp(2) = X_schematic(2,j,k) + du/4; 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

        x_tmp(1) = X_schematic(1,j,k) + 3*du/4; 
        x_tmp(2) = X_schematic(1,j,k) + du; 
        y_tmp(1) = X_schematic(2,j,k) + du/4; 
        y_tmp(2) = X_schematic(2,j,k); 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

    end 
end 

j=1; 
for k=1:k_max
    if is_internal(j,k) 
        x_tmp(1) = X_schematic(1,j,k); 
        x_tmp(2) = X_schematic(1,j,k) - du/2; 
        y_tmp(1) = X_schematic(2,j,k); 
        y_tmp(2) = X_schematic(2,j,k); 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

        % make a little zig zag to show periodicity 
        x_tmp(1) = X_schematic(1,j,k) - du/2; 
        x_tmp(2) = X_schematic(1,j,k) - 3*du/4; 
        y_tmp(1) = X_schematic(2,j,k); 
        y_tmp(2) = X_schematic(2,j,k) - du/4; 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

        x_tmp(1) = X_schematic(1,j,k) - 3*du/4; 
        x_tmp(2) = X_schematic(1,j,k) - du; 
        y_tmp(1) = X_schematic(2,j,k) - du/4; 
        y_tmp(2) = X_schematic(2,j,k); 
        plot(x_tmp, y_tmp, plot_options, 'MarkerSize', marker_size); 

    end 
end 




axis equal;    

x_min = -2*du; 
x_max = 1 + 2*du; 
y_min = .7; 
y_max = 1 + du; % add one du here for a little space 


axis([x_min x_max y_min y_max])

set(gcf,'color',[1 1 1])
xlabel('u')
ylabel('v')

printfig(fig, 'mesh_schematic_leaflet'); 








