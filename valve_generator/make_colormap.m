function cmap = make_colormap(len, extended)
% 
% returns a custom color map of length len 
% 
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

try_number = 0; 

if try_number == 1 

    hue = linspace(3/6,0,len)'; 

    min_full_saturation = floor(len * 3/4); 

    saturation = linspace(0,1,min_full_saturation)'; 
    saturation = [saturation; ones(len-min_full_saturation,1);]; 

    max_full_value = floor(len * 3/4); 
    value = ones(max_full_value, 1); 
    value = [value; linspace(1,0,len-max_full_value)']; 

    hsv = [hue, saturation, value]; 

    cmap = hsv2rgb(hsv); 

elseif try_number == 2
    
    min_inverted_hot = floor(len/2);
    
    hue        = linspace(4/6,1/6,min_inverted_hot)'; 
    saturation = linspace(0,1,min_inverted_hot)'; 
    value      = ones(min_inverted_hot, 1); 
    
    hsv = [hue, saturation, value]; 
    
    
    hot_values = hot(len - min_inverted_hot);
    hot_values = flipud(hot_values);
    hot_hsv = rgb2hsv(hot_values); 
    
    hsv = [hsv; hot_hsv]; 
    
    cmap = hsv2rgb(hsv);
    
else 
    
    % hue goes linearly from some value
    % to zero, then is zero (red) near the end 
    % 
    % saturation goes linearly to one then stays there 
    % (more grays as you go up)
    % 
    % value is one, then linearly to zero 
    % bright at first, then darker 
    
    if exist('extended', 'var') && extended 
        
        % put some dark colors at the top to ignore
        % linear things go up earlier to get reds closer to middle
        
        linear_hue_frac = .65; 
        max_hue = 4/6; 
    
        linear_saturation_frac = .6; 
    
        value_flat_frac = .5; % .84;     

    else 
    
        linear_hue_frac = .85; 
        max_hue = 4/6; 

        linear_saturation_frac = .6; 

        value_flat_frac = .84;     
    end 
    
    sloped_hue_length = floor(linear_hue_frac * len); 
    hue = linspace(max_hue,0,sloped_hue_length)'; 
    hue = [hue; zeros(len - sloped_hue_length,1)]; 

    sloped_saturation_length = floor(linear_saturation_frac * len); 
    saturation = linspace(0,1,sloped_saturation_length)'; 
    saturation = [saturation; ones(len - sloped_saturation_length,1)]; 
    
    flat_value_length = floor(value_flat_frac * len); 
    value = ones(flat_value_length, 1); 
    value = [value; linspace(1,0,len - flat_value_length)']; 
    
    hsv = [hue, saturation, value]; 
    
    cmap = hsv2rgb(hsv);
end 


debug = false; 
if debug 
    figure; 
    plot(hsv(:,1)); 
    title('hue')

    figure; 
    plot(hsv(:,2)); 
    title('saturation')

    figure; 
    plot(hsv(:,3)); 
    title('value')
end 

    
    