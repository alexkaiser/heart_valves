function [] = output_valve_ring(valve)

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

X     = valve.leaflets(1).X; 
k_max = valve.leaflets.k_max; 

x_ring = squeeze(X(1,:,k_max)); 
y_ring = squeeze(X(2,:,k_max)); 

fig = figure; 
plot(x_ring, y_ring, 'k')
hold on 

% patch periodic
x_patch = [x_ring(end), x_ring(1)]; 
y_patch = [y_ring(end), y_ring(1)]; 

plot(x_patch, y_patch, 'k'); 


% add two little ticks 
tan = [x_ring(end) - x_ring(2), y_ring(end) - y_ring(2)]; 

normal = [tan(2) -tan(1)]; 
normal = normal / norm(normal); 

len = .15; 

x_tick = [x_ring(1) + len*normal(1), x_ring(1) - len*normal(1)]; 
y_tick = [y_ring(1) + len*normal(2), y_ring(1) - len*normal(2)]; 
plot(x_tick, y_tick, 'k') 


idx_right = valve.leaflets(1).j_range_anterior(end); 
tan = [x_ring(idx_right+1) - x_ring(idx_right-1), y_ring(idx_right+1) - y_ring(idx_right-1)]; 

normal = [tan(2) -tan(1)]; 
normal = normal / norm(normal); 

x_tick = [x_ring(idx_right) + len*normal(1), x_ring(idx_right) - len*normal(1)]; 
y_tick = [y_ring(idx_right) + len*normal(2), y_ring(idx_right) - len*normal(2)]; 
plot(x_tick, y_tick, 'k') 





buffer = .3; 

axis equal; 
axis([min(x_ring)-buffer max(x_ring)+buffer min(y_ring)-buffer max(y_ring)+buffer])

xlabel('cm')
ylabel('cm')
title('Valve ring'); 

set(gcf,'color',[1 1 1])
printfig(fig, 'valve_ring_2d.eps'); 

