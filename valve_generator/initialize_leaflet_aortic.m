function [leaflet valve] = initialize_leaflet_aortic(name,                         ...
                                                     N,                            ...
                                                     tension_coeffs,               ...
                                                     p_0,                          ...
                                                     valve)
%
% Builds leaflet data structures 
% 
% Input:                                   
%                                   
%     N                             Size parameter for leaflet 
% 
% Output 
% 
%     leaflet                       Fully initialized leaflet data structure 
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
  
leaflet.name               = name; 
leaflet.N                  = N; 

leaflet.decreasing_tension = valve.decreasing_tension;

leaflet.diff_eqns = valve.diff_eqns; 
leaflet.jacobian  = valve.jacobian;

leaflet.skeleton = valve.skeleton; 


if isfield(valve, 'careful_early_steps')
    leaflet.careful_early_steps = valve.careful_early_steps; 
    if leaflet.careful_early_steps
        
        if isfield(valve, 'careful_early_step_coeff')
            leaflet.careful_early_step_coeff = valve.careful_early_step_coeff; 
        else
            warning('Using default careful_early_step_coeff of 1/2')
            leaflet.careful_early_step_coeff = 1/2; 
        end
        
        if isfield(valve, 'residual_decrease_to_double')
            leaflet.residual_decrease_to_double = valve.residual_decrease_to_double; 
        else 
            warning('Using default residual_decrease_to_double of 1/2')
            leaflet.residual_decrease_to_double = 1/2; 
        end 
    end     
end 

if isfield(valve, 'variety')
    leaflet.variety = valve.variety; 
end 

if isfield(leaflet, 'variety') && strcmp(leaflet.variety, 'bicuspid') && (mod(N,2) ~= 0)
    error('Bicuspid aortic valve N must be a multiple of two'); 
end
    
if mod(N,3) ~= 0
    error('Aortic valve N must be a multiple of three'); 
end 


leaflet.du = 1/N; 

leaflet.r = valve.r; 

leaflet.ds = valve.ds; 

% information about geometry 
leaflet = get_util_arrays_aortic(leaflet, valve); 

% build actual data structure 
leaflet.X = build_initial_fibers_aortic(leaflet, valve); 

% Scalar pressure to support 
leaflet.p_0 = p_0; 

% Total number of internal leaflet coordinates (three times number of vertices)
leaflet.total_internal_leaflet    = 3*sum(leaflet.is_internal(:)); 
    
% set coefficients on tensions
[leaflet, valve] = set_tension_coeffs(leaflet, valve, tension_coeffs); 

% parameter structure for collagen based nonlinear constitutive 
leaflet.collagen_curve             = get_collagen_curve_parameters(); 
leaflet.collagen_constitutive_circ = valve.collagen_constitutive_circ; 
leaflet.collagen_constitutive_rad  = valve.collagen_constitutive_rad; 



