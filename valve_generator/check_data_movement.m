
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

N = 4; 

% Main data structure with everything 
valve.N = N; 
valve.tol_global = 1e-10;
valve.max_it = 40; 


% Valve skeleton parameters 
valve.r = 1.606587877768772; 
valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 
valve.split_papillary = false; 
valve.radial_and_circumferential = true; 

% function pointers 
valve.diff_eqns = @difference_equations; 
valve.jacobian  = @build_jacobian; 

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

% box width 
valve.L = 2.5; 

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
valve.pressure_tension_ratio = 1.5; 


% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
valve.refinement = N/32.0; 

valve.p_physical = 100; 

% scaling for target points 
valve.target_multiplier = 40; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
valve.n_lagrangian_tracers = 8; 

% Uses configuration of X 
valve.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
valve.num_copies = 3; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_springs_leaflet = false; 


% posterior leaflet data structure 
reflect_x   = true; 
total_angle = pi + pi/6 + pi/12; 

% filter paramters
a = 1.0; 
h = 2.0; 
r = valve.r; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
radial_and_circumferential = true; 

% This many points are placed in a flat manner 
% around the point of the free edge 
trapezoidal_flat_points = 0; % N/8; 

% Spring constants in two directions 
alpha    = 1.0; 
beta     = 1.0; 
p_0      = 0.0; % no pressure for now 
ref_frac = 0.7; % generic spring constants reduced by this much 

% Chordae parameters 
k_0          = 2.0; 
k_multiplier = 2.0; 
tree_frac    = 0.5;



valve.posterior = initialize_leaflet(N,                            ... 
                                     reflect_x,                    ... 
                                     total_angle,                  ...    
                                     a,                            ... 
                                     h,                            ... 
                                     r,                            ... 
                                     valve.left_papillary,         ... 
                                     valve.right_papillary,        ... 
                                     radial_and_circumferential,   ...  
                                     alpha,                        ... 
                                     beta,                         ... 
                                     p_0,                          ... 
                                     ref_frac,                     ...  
                                     k_0,                          ... 
                                     k_multiplier,                 ... 
                                     tree_frac,                    ...
                                     trapezoidal_flat_points); 

'version with no movement, including b.c.s'
valve.posterior.X 


% 'linear order on internal points'
X_linearized = linearize_internal_points(valve.posterior, valve.posterior.X); 


'moved back to 3d ordering'
leaflet = internal_points_to_2d(X_linearized, valve.posterior);  
leaflet.X 


'with chordae before linearization'
valve.posterior.X 
valve.posterior.chordae.C_left
valve.posterior.chordae.C_right 

'linearized with chordae'
X_and_chordae_linearized = linearize_internal_points(valve.posterior, valve.posterior.X, valve.posterior.chordae.C_left, valve.posterior.chordae.C_right) 


'after return to normal data structure'
valve.posterior = internal_points_to_2d(X_and_chordae_linearized, valve.posterior);  
valve.posterior.X 
valve.posterior.chordae.C_left
valve.posterior.chordae.C_right 











