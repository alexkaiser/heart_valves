function [valve] = initialize_valve_data_structures_hocm_d(N, attached, leaflet_only, optimization, decreasing_tension)
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
valve.careful_early_steps         = false; 
valve.careful_early_step_coeff    = 1/8; 
valve.residual_decrease_to_double = 1/2; 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 0;  
valve.max_total_fails       = 0; 

if exist('attached', 'var') 
    valve.attached = attached; 
else 
    valve.attached = false; 
end 

split_papillary = true; 
valve.split_papillary = split_papillary; 
valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.leaflet_only = leaflet_only; 
valve.optimization = optimization; 

valve.decreasing_tension = decreasing_tension; 


valve.diff_eqns = @difference_equations_bead_slip; 
valve.jacobian  = @build_jacobian_bead_slip;

valve.targets_for_bcs = false; 
valve.targets_for_bcs_ref_only = false; 

% general solve parameters

% does not place partition
% mitral ring and papillary are still target points 
valve.in_heart = true; 

% name 
if valve.in_heart
    valve.base_name = sprintf('mitral_no_partition_%d', N); 
    valve.extra_radius_hoops = 0; 0.13; % adds points out the partition up to this amount 
                                     % two extra hoops at 256
else 
    valve.base_name = sprintf('mitral_tree_%d', N); 
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
valve.num_copies = 3; 

% add flags to spring files 
% to view and output with a stride 


valve.output.leaflets       = [1;1;1]; 
valve.output.stride_leaflet = max(1,N/128); 
valve.output.chordae        = [1;1;1]; 
valve.output.mesh           = [1;0;0]; 
valve.output.cartesian_mesh = [0;0;0]; 
valve.output.stride_mesh    = N/32; 

valve.papillary_movement_times = [0 .1 .46 .5 .8]; 


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

% name of structure 
name = 'leaflet'; 

% commissural tree version 
% but without explicit commissural leaflets 

valve.p_physical = 120 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

% valve.dip_anterior_systole = true; 
% valve.r_dip = 0.75; 
% valve.total_angle_dip = pi; 

valve.L = 3.0; 

valve.skeleton = get_skeleton_hcm_d(); 
valve.r = valve.skeleton.r; 

valve.diastolic_increment = [0.0; 0.0; 0.0]; 


% Base constants, individual pieces are tuned relative to these values

% tension coefficients structure 

% pressure / tension coefficient ratio
% this tension coefficient is the maximum tension that a fiber can support
% valve.pressure_tension_ratio = 0.055; % 0.11 * 0.975; 
tension_coeffs.pressure_tension_ratio = 0.041; 

tension_coeffs.dec_tension_coeff_base = 7; 


% max tensions in leaflets 
tension_coeffs.alpha_anterior       = 1.75;  % circumferential 
tension_coeffs.beta_anterior        = 1.2;  % radial
tension_coeffs.alpha_posterior      = 1.0;  % circumferential 
tension_coeffs.beta_posterior       = 0.75;  % radial
tension_coeffs.alpha_hoops          = 0.5;  % circumferential hoops     
tension_coeffs.alpha_edge_connector = 1.75;  % circumferential free edge connector 
tension_coeffs.beta_edge_connector  = 0.01;  % circumferential free edge connector


% decreasing tension coefficients 
tension_coeffs.c_circ_dec_anterior       = 5.0;  % circumferential 
tension_coeffs.c_rad_dec_anterior        = 1.5;  % radial
tension_coeffs.c_circ_dec_posterior      = 3.0;  % circumferential 
tension_coeffs.c_rad_dec_posterior       = 0.47;  % radial
tension_coeffs.c_circ_dec_hoops          = 2.0;  % circumferential hoops
tension_coeffs.c_rad_dec_hoops_anterior  = 0.5;  % radial hoops, anterior part 
tension_coeffs.c_rad_dec_hoops_posterior = 0.3;  % radial hoops, posterior part 
tension_coeffs.c_dec_tension_chordae     = 1.0;  % chordae

tension_coeffs.c_circ_dec_edge_connector  = 1.0;  % circumferential hoops
tension_coeffs.c_rad_dec_edge_connector   = 2.0;  % circumferential hoops

% places this many periodic rings above 
n_rings_periodic = max(1,N/64); 

% places circumferential fibers this many below hoops 
% if the location is not already covered by leaflet
n_edge_connectors = max(1,N/64);  

% No explicit commissural leaflet here 
N_anterior = N/2; 

angles.anterior = 5*pi/6; 

% Posterior takes whatever is left 
N_posterior = N - N_anterior; 

% store these 
valve.N_anterior   = N_anterior; 
valve.N_posterior  = N_posterior;
valve.commissural_leaflets = false; 
valve.N_commissure = 0; 


N_per_direction   = [N_anterior/2, N_anterior/2, N_posterior/2, N_posterior/2]; 

% Anterior goes down then up 
leaflet_direction = [-1, 1]; 

% Posterior goes down then up 
leaflet_direction = [leaflet_direction, -1, 1]; 

% No offset, starting at commissure 
leaflet_N_start = 0; 

% changes entire tree strength by constants 
tension_coeffs.tree_tension_multiplier = 0.75; 

% Leaf tensions are all modified 
tension_coeffs.leaf_tension_base = .9; 

% Base total root tension 
% The value 0.5905 works well on each tree when using separate solves and two leaflets 
% Controls constant tension at the root of the tree 
tension_coeffs.root_tension_base = .9 * 0.5905; 

% this array determines the fraction of N_orig which each tree takes up 
% this allows us to determine initial fractions of constants that go to each tree 
frac_of_n_orig = [1/ 4; 1/ 4;  ...   % anterior  
                  1/16; 1/16;  ...   % comm
                  1/ 8; 1/ 8;  ...   % posterior
                  1/16; 1/16];       % comm

% change these to manipulate individial tree coefficients 
% for sanity reasons, these shuold mostly be one unless you have a good reason to change 
% note that these are scaled by the fraction of the leaflet that they take up 
k_0_1_coeff    = frac_of_n_orig .*    ... 
                 [2.2; 2.2;           ...       % anterior  
                  2.0; 2.0;           ...       % anterior and comm, comm and posterior       
                  1.6; 1.6;           ...       % posterior
                  2.0; 2.0];                    % posterior and comm, comm and anterior

k_root_coeff   = frac_of_n_orig .*    ... 
                [ 2.2; 2.2;           ...       % anterior  
                  2.0; 2.0;           ...       % anterior and comm, comm and posterior       
                  1.6; 1.6;           ...       % posterior
                  2.0; 2.0];                    % posterior and comm, comm and anterior


tension_coeffs.c_dec_chordae_leaf = (1/N)  * [1.0; 1.0; ...       % anterior  
                                              1.0; 1.0; ...       % comm
                                              1.0; 1.0; ...       % posterior
                                              1.0; 1.0];          % comm 

% root constants do not scale, because the root 
% should maintain a consistent length when mesh is changed 
tension_coeffs.c_dec_chordae_root = [0.05;  0.05; ...       % anterior  
                                     1/256; 1/256; ...       % comm
                                     1/256; 1/256; ...       % posterior
                                     1/256; 1/256];          % comm 

n_leaves = N * frac_of_n_orig; 

% concatenate all relevant arrays
% n_leaves           = [n_leaves_anterior; n_leaves_posterior_and_comm];
% papillary          = [papillary_anterior, papillary_posterior_and_comm]; 
% k_0_1              = [k_0_1_anterior; k_0_1_posterior_and_comm]; 
% k_root             = [k_root_anterior; k_root_posterior_and_comm]; 

tension_coeffs.k_0_1  = k_0_1_coeff; 
tension_coeffs.k_root = k_root_coeff; 
    
% default value for tree offset 
if ~exist('tree_n_start', 'var')
    tree_n_start = 1; 
end 


% valve.r        = valve.skeleton.r; 

% scaling for target points 
% note that this does not include copies 
% and scaling for copies is handled by the output routine 

% scales for by mesh width for consistant total mesh force on ring 
valve.target_net_unscaled       = 8 / valve.N; 

% does not scale since total number of points is constant 
valve.target_papillary_unscaled = 10 * 40/128; 

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


[leaflet valve] = initialize_leaflet_bead_slip(name,                                ... 
                                               N,                                   ...
                                               reflect_x,                           ... 
                                               angles,                              ...
                                               valve.skeleton.papillary,            ... 
                                               n_leaves,                            ...
                                               tree_n_start,                        ...
                                               leaflet_direction,                   ...
                                               leaflet_N_start,                     ...
                                               N_per_direction,                     ...
                                               radial_and_circumferential,          ...  
                                               tension_coeffs,                      ... 
                                               p_0,                                 ... 
                                               n_rings_periodic,                    ...
                                               n_edge_connectors,                   ... 
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



