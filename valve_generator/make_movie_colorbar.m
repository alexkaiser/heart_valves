function make_movie_colorbar(max_velocity, format_string, time, file_name)

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

% process speeds to colormap 
all_colors = true; 
n_colors = 500;
colormap_temp = make_colormap(n_colors); 
if all_colors 
    range = 1:floor(   length(colormap_temp)); 
else 
    range = 1:floor(.9*length(colormap_temp)); 
end 
colormap_croppeed = colormap_temp(range,:); 


% colorbar for movie, bigger and bolder 
n_ticks = 4; 
tick_array = linspace(0,1,n_ticks); 
tick_labels = {}; 
for i=1:length(tick_array)
    tick=tick_array(i); 
    v = tick * max_velocity; 
    tick_labels{i} = sprintf('%d', v); 
end 

fontsize = 24 * 4;

fig = figure('visible','off'); 
set(fig, 'Renderer', 'Painters');
set(fig, 'Position', [0 0 444 4320])
set(fig, 'PaperPositionMode','auto')

colormap(colormap_croppeed); 
cbar = colorbar('Ticks', tick_array, 'TickLabels', tick_labels, 'fontweight', 'bold', 'fontsize',fontsize); 

if strcmp(format_string, '.eps')
    cbar.Label.String = sprintf('\n\nt (s)\n%.2f\n\n\n|u|\n(cm/s)', time);
else
    cbar.Label.String = sprintf('t (s)\n%.2f\n\n\n|u|\n(cm/s)', time);
end 

cbar.Label.FontSize = fontsize; 
cbar.Label.Position = [0 1.05]; 
cbar.Label.Rotation = 0;
cbar.Label.FontWeight = 'bold'; 
cbar.TickDirection = 'out'; 
cbar.AxisLocation = 'in'; 


% position and sizing of the bar itself relative to the block 
% does not include labels
cbar.Position = [0.6 0.27 .4 .5]; 


grid off 
axis off 

print(fig, format_string, file_name, '-r0'); 


