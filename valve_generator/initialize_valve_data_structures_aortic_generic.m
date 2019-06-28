function [valve] = initialize_valve_data_structures_aortic_generic(N)
% 
% Initializes data structures for full solve.  
% 
% Parameters are declared here.
% Should be a script, but want to return in the structures 
% 
% Input: 
%     N   Size parameter used throughout 
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

% Main data structure with everything 
valve.N = N; 

% effective infinity by default 
valve.max_it                = 1e8; 
valve.max_continuations     = 1e8; 

% shirinks initial 
valve.careful_early_steps         = true; 
valve.careful_early_step_coeff    = 1/8; 
valve.residual_decrease_to_double = 1/2; 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 0;  
valve.max_total_fails       = 0; 


valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.decreasing_tension = true; 

warning('fix diff eqns for aortic')
valve.diff_eqns = @difference_equations_aortic; 
valve.jacobian  = @build_jacobian_aortic;

valve.targets_for_bcs = false; 
valve.targets_for_bcs_ref_only = false; 

% general solve parameters
name = 'aortic'; 
valve.name = name; 

% does not place partition
valve.in_heart = true; 

% name 
if valve.in_heart
    valve.base_name = sprintf('aortic_no_partition_%d', N); 
    valve.extra_radius_hoops = 0.13; % adds points out the partition up to this amount 
                                     % two extra hoops at 256
else 
    valve.base_name = sprintf('aortic_%d', N); 
end 
MMHG_TO_CGS     = 1333.22368;


% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
valve.n_lagrangian_tracers = 8; 

% Uses configuration of X 
valve.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
valve.num_copies = 1; 

% add flags to spring files 
% to view and output with a stride 

valve.output.leaflets       = [1;1;1]; 
valve.output.stride_leaflet = max(1,N/128); 
valve.output.mesh           = [1;0;0]; 
valve.output.cartesian_mesh = [0;0;0]; 
valve.output.stride_mesh    = N/32; 


% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_constitutive = true; 

% Constant strain of pressurized configuration 
valve.strain = .16; 

% no reflections in this version 
reflect_x = false; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
% Always true in this version `
radial_and_circumferential = true; 

% physical units create a scalar multiple of the old 
% this multiple is large number, so we want to scale the old tolerance accordingly 
% 8.3326e-04 is a good number here
valve.tol_global = 1e-3;


% commissural tree version 
% but without explicit commissural leaflets 
valve.p_physical = 40 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

% valve.dip_anterior_systole = true; 
% valve.r_dip = 0.75; 
% valve.total_angle_dip = pi; 

valve.L = 3.0; 

valve.skeleton = get_skeleton_aortic_generic(); 
valve.r = valve.skeleton.r; 

valve.diastolic_increment = [0.0; 0.0; 0.0]; 


% Base constants, individual pieces are tuned relative to these values

% tension coefficients structure 

% pressure / tension coefficient ratio
% this tension coefficient is the maximum tension that a fiber can support
% valve.pressure_tension_ratio = 0.055; % 0.11 * 0.975; 
tension_coeffs.pressure_tension_ratio = 0.055; 

tension_coeffs.dec_tension_coeff_base = 4.6; 


% max tensions in leaflets 
tension_coeffs.alpha = 1.0;   % circumferential 
tension_coeffs.beta  = 0.05;  % radial

% decreasing tension coefficients 
tension_coeffs.c_circ_dec       = 1.0;  % circumferential 
tension_coeffs.c_rad_dec        = 1.0;  % radial

% scaling for target points 
% note that this does not include copies 
% and scaling for copies is handled by the output routine 

% scales for by mesh width for consistant total mesh force on ring 
valve.target_net_unscaled       = 8 / valve.N; 

% does not scale since total number of points is constant 
valve.target_papillary_unscaled = 2 * 40/128; 

% viscoelastic damping coefficients for net, does not include copies 
valve.eta_net_unscaled = 0.0; %valve.target_net_unscaled/5000; 

% viscoelastic damping coefficients for root attachments, does not include copies  
valve.eta_papillary_unscaled = 0.0; valve.target_papillary_unscaled/500; 

% if nonzero, linear springs of rest length with spacing between the layers 
% are placed with this value 
% final formula is multiplied by valve.tension_base  
valve.kappa_cross_layer_multipler = 1e4 / 256^2; 

% Approximate Lagrangian mesh spacing at ring 
% Used for later splitting of springs 
% If any spring is placed at more than double this length an extra vertex is placed
valve.ds = 2*pi*valve.skeleton.r / N; 

[leaflet valve] = initialize_leaflet_aortic(name,                                ... 
                                            N,                                   ...
                                            tension_coeffs,                      ... 
                                            p_0,                                 ... 
                                            valve);  

valve.leaflets(1) = leaflet; 
    

% viscoelastic damping coefficients springs 
% eta, damping coeff here, is multiplied by the coefficient on the 
% associated spring 
% note that linear springs and collagen springs have vastly different constants 
% and these are tuned manually to make the dashpot constants equal order of magnitude
valve.eta_multiplier_linear   = 0; 
valve.eta_multiplier_collagen = 0; 

valve_plot(valve); 
pause(.1); 

disp('Done with initialize.'); 

