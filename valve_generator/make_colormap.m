function cmap = make_colormap(len, extended, export_xml)
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

if ~exist('export_xml', 'var')
    export_xml = false; 
end 

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

if export_xml 
    f = fopen('colormap.xml', 'w');
    
    n_colors = size(cmap,1); 
    
    x_vals = linspace(0,1,n_colors);

    fprintf(f, "<ColorMaps>\n");
    fprintf(f, '<ColorMap name="RainbowLinear" space="RGB" interpolationspace="RGB" interpolationtype="linear" creator="CCC-Tool">\n');
    
    for i=1:n_colors
        fprintf(f, '<Point x="%.14f" o="1" r="%.14f" g="%.14f" b="%.14f" cms="1" isMoT="true"/>\n', x_vals(i), cmap(i,1), cmap(i,2), cmap(i,3));
    end 
    
    fprintf(f, '<NaN r="0" g="0" b="0"/>\n');
    fprintf(f, '<Above r="0" g="0" b="0"/>\n');
    fprintf(f, '<Below r="1" g="1" b="1"/>\n');
    fprintf(f, '</ColorMap>\n');
    fprintf(f, '</ColorMaps>\n');
    
    
%     <ColorMaps>
%     <ColorMap name="YellowGreenPurpleBrown" space="RGB" interpolationspace="lab" interpolationtype="linear" creator="CCC-Tool">
%     <Point x="0" o="1" r="0.741" g="0.9116833333333333" b="0.95" cms="1" isMoT="true"/>
%     <Point x="0.125" o="1" r="0.5192000000000001" g="0.7236533333333334" b="0.88" cms="1" isMoT="true"/>
%     <Point x="0.42" o="1" r="0.246" g="0.3639999999999999" b="0.6" cms="1" isMoT="true"/>
%     <Point x="0.5" o="1" r="0.06660000000000002" g="0.11210999999999985" b="0.37" cms="1" isMoT="true"/>
%     <Point x="0.5" o="1" r="0.17359166666666662" g="0.37" b="0.03329999999999999" cms="1" isMoT="true"/>
%     <Point x="0.65" o="1" r="0.2943791557480731" g="0.59" b="0.21239999999999998" cms="1" isMoT="true"/>
%     <Point x="1" o="1" r="0.8486949532130028" g="1" b="0.71" cms="1" isMoT="true"/>
%     <NaN r="0" g="0" b="0"/>
%     <Above r="0" g="0" b="0"/>
%     <Below r="0" g="0" b="0"/>
%     </ColorMap>
%     </ColorMaps>
end 





    
    