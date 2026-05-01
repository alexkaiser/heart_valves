function [valve] = init_aortic_historical_3(N)
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

valve.rotate_identical_leaflets = false;

% effective infinity by default 
valve.max_it                = 1e8; 
valve.max_continuations     = 1e8; 

% shirinks initial 
valve.careful_early_steps         = false; 
if valve.careful_early_steps
    valve.careful_early_step_coeff    = 1/8; 
    valve.residual_decrease_to_double = 1/2; 
end 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 0;  
valve.max_total_fails       = 0; 


valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.decreasing_tension = true; 

valve.diff_eqns = @difference_equations_aortic; 
valve.jacobian  = @build_jacobian_aortic;

% general solve parameters
name = 'aortic'; 
valve.name = name; 

variety= 'bicuspid'; 
valve.variety = variety; 

% does not place partition
valve.in_heart = true; 

valve.inst_file_ring = false; 

valve.base_name = sprintf('aortic_no_partition_%d', N); 
valve.extra_radius_hoops = 0.0; % adds points out the partition up to this amount 


valve.tight_cylinder = true; 
valve.z_extra_cylinder = 0.3; 

% for normal_3
% if starting at origin 
% th = 2*pi/3; 
% valve.initial_translation_aortic = 0.1 * [cos(th); sin(th); 0]; 
% valve.initial_translation_aortic = valve.initial_translation_aortic  + 0.1 * [1; 0; 0]; 
% valve.initial_rotation_aortic = rotation_matrix_z(pi/3 + pi/6 + pi/12);
% valve.transformation_vertex_file = 'historical_3_vbr.vertex';



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

% valve.copy_spring_weights = [1/2 1/4 1/4];

valve.normal_thicken = true; 
% nominal aortic valve thickness
valve.normal_thickness = 0.044; 

valve.extrusion_out = true;

% respace on annulus in 3d 
% if false, spaced wrt theta 
valve.annulus_points_even_spacing = true; 

valve.annulus_to_comm = true; 


valve.dirichlet_free_edge = false; 

valve.dirichlet_free_edge_with_ref_only = true; 




% provides a bending resistance for the final solve 
% for initial conditions 
% this is just to get a reasonable initial condition 
valve.k_bend_radial_ref_only = 0; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_constitutive_circ = 'aortic_circ'; 
valve.collagen_constitutive_rad  = 'aortic_rad'; 
 
% Constant strain of pressurized configuration 
valve.strain_circ = .15; 
valve.strain_rad  = .54; 

valve.extra_stretch_radial_dirichlet_free_edge = 1.0 * valve.strain_rad + 1.0; 

% physical units create a scalar multiple of the old 
% this multiple is large number, so we want to scale the old tolerance accordingly 
% 8.3326e-04 is a good number here
valve.tol_global = 1e-3;


% commissural tree version 
% but without explicit commissural leaflets 
valve.p_physical = 60 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

valve.p_final = 0.0 * MMHG_TO_CGS;  


valve.dirichlet_free_edge_comm_ref_only = false; 
valve.n_fixed_comm = max(1, floor(8*N/192));  
valve.p_final_fixed_comm = 0.1 * MMHG_TO_CGS;  


valve.L = 2.25; 

% r_stj = 2.5/2; % 25 mm valve 
% r_temp = r_stj; 
% hc = 0.5 * r_stj; 
% h1 = 1.4 * r_stj - hc; 
% r_commissure = r_stj; 
% % place the post only if not using the full annulus geometry 
% place_vertical_post = ~valve.annulus_to_comm;
% valve.skeleton = get_skeleton_aortic_generic(r_temp, h1, hc, r_commissure, place_vertical_post); 
% % valve.skeleton = get_skeleton_aortic_generic(); 
% valve.r = valve.skeleton.r; 

valve.skeleton = get_skeleton_aortic_hist3();

valve.r = valve.skeleton.r; 


% valve_ring_pts_raw_LR_cusp = 0.1 * [ 
% -16.983776092529297, -161.2555389404297, -163.5882568359375,
% -13.077391624450684, -157.9479522705078, -169.11514282226563,
% -7.8202805519104, -158.03750610351563, -173.54421997070313,
% -2.486410140991211, -160.36691284179688, -175.14462280273438,
% 1.3799179792404175, -165.36549377441406, -175.52175903320313,
% -1.3226488828659058, -174.33387756347656, -178.4099578857422,
% -9.562684544740621, -181.61004635596518, -180.93456165010247,
% -15.33113956451416, -183.05392456054688, -179.898681640625,
% -19.65570831298828, -183.8021697998047, -177.125,
% -23.952800750732422, -184.1968994140625, -170.5049285888672]';
% 
% % these happen to be in right to left order from segmentation 
% valve_ring_pts_raw_LR_cusp = fliplr(valve_ring_pts_raw_LR_cusp);

% consider aligning before model construction 
% R_0 = valve.initial_rotation_aortic; 
% T_0 = valve.initial_translation_aortic; 
% apply_inverse = true;
% valve.skeleton.valve_ring_pts = coordinate_transformation_vertices(valve_ring_pts_raw_LR_cusp, valve.transformation_vertex_file, R_0, T_0, apply_inverse);
% valve.skeleton.valve_ring_pts = valve_ring_pts_raw_LR_cusp;   
% plot3(valve.skeleton.valve_ring_pts(1,:), valve.skeleton.valve_ring_pts(2,:), valve.skeleton.valve_ring_pts(3,:),'*-');


valve.place_cylinder = false; 
valve.z_max_cylinder = (pi/3) * valve.r; 
valve.z_min_cylinder = 0.0; 


valve.n_layers_cylinder = 3; 


% comm_raise_normal_height = 0.8 * valve.skeleton.r * 2; 


% Base constants, individual pieces are tuned relative to these values

% tension coefficients structure 

% LR leaflet 
tension_coeffs_lr.pressure_tension_ratio = 0.00477; 
tension_coeffs_lr.dec_tension_coeff_base = 20.0; 

% max tensions in leaflets 
tension_coeffs_lr.alpha = 1.6;   % circumferential 
tension_coeffs_lr.beta  = 0.055;   % radial

% decreasing tension coefficients 
tension_coeffs_lr.c_circ_dec       = 3.15;  % circumferential 
tension_coeffs_lr.c_rad_dec        = 0.9;  % radial
tension_coeffs_lr.c_circ_dec_annulus = 1.91;
tension_coeffs_lr.c_circ_dec_free_edge_percentage = 0.0;

% non coronary leaflet
tension_coeffs_non.pressure_tension_ratio = 0.00477; 
tension_coeffs_non.dec_tension_coeff_base = 20.0; 

% max tensions in leaflets 
tension_coeffs_non.alpha = 1.6;   % circumferential 
tension_coeffs_non.beta  = 0.055;   % radial

% decreasing tension coefficients 
tension_coeffs_non.c_circ_dec       = 3.22;  % circumferential 
tension_coeffs_non.c_rad_dec        = 0.84;  % radial
tension_coeffs_non.c_circ_dec_annulus = 1.91;
tension_coeffs_non.c_circ_dec_free_edge_percentage = 0.0;


% scaling for target points 
% note that this does not include copies 
% and scaling for copies is handled by the output routine 

% scales for by mesh width for consistant total mesh force on ring 
valve.target_net_unscaled       = (8 / valve.N) * (192/N); 

% does not scale since total number of points is constant 
valve.target_papillary_unscaled = 2 * 40/128; 

% viscoelastic damping coefficients for net, does not include copies 
valve.eta_net_unscaled = 0; % 1e-5 * valve.target_net_unscaled; 

% viscoelastic damping coefficients for root attachments, does not include copies  
valve.eta_papillary_unscaled = 0.0; 

% if nonzero, linear springs of rest length with spacing between the layers 
% are placed with this value 
% final formula is multiplied by valve.tension_base  
valve.kappa_cross_layer_multipler = 2 * (384/N)^2 * 1e4 / 256^2;

% valve.k_bend_radial = [0 0 1e5 1e5] * 192/N;
valve.k_bend_radial = 1e4 * 192/N;
% valve.k_bend_radial_annulus = 1e2 * 192/N;
valve.k_bend_radial_free_edge = 0; 
valve.k_bend_radial_free_edge_percentage = 0; 
valve.k_bend_circ = 1e4 * 192/N;
valve.k_bend_circ_free_edge = 0; 
valve.k_bend_circ_free_edge_percentage = 0;

valve.k_bend_cross_layer = 1e4 * 192/N;



% coaptation height for Swanson and Clark 
% width is the same 
% valve.k_bend_nodule_length = 0.17*2*valve.skeleton.r;
% valve.k_bend_nodule        = 1e5 * 192/N;

% valve.k_bend_radial_interp_pts  = [0    .36  .46   1];

% valve.kappa_radial_free_edge_compressive_unscaled = 1e3 / 256^2;
% valve.kappa_radial_free_edge_compressive_percentage = 0.4;
% valve.kappa_radial_free_edge_compressive_stretch = 1.54;
% valve.kappa_radial_free_edge_compressive_fn_idx = 4;

% Approximate Lagrangian mesh spacing at ring 
% Used for later splitting of springs 
% If any spring is placed at more than double this length an extra vertex is placed
valve.ds = 2*pi*valve.skeleton.r / N; 

[leaflet_lr valve] = initialize_leaflet_aortic(name,                             ... 
                                            N,                                   ...
                                            tension_coeffs_lr,                   ... 
                                            p_0,                                 ... 
                                            valve,                               ...
                                            valve.skeleton.ring_pts_LR_cusp_model_coords);  

valve.leaflets(1) = leaflet_lr;


[leaflet_non valve] = initialize_leaflet_aortic(name,                            ... 
                                            N,                                   ...
                                            tension_coeffs_non,                  ... 
                                            p_0,                                 ... 
                                            valve,                               ...
                                            valve.skeleton.ring_pts_Non_cusp_model_coords);  

valve.leaflets(2) = leaflet_non;



if isfield(valve, 'extrusion_out') 
    for leaflet_idx = 1:length(valve.leaflets)
        valve.leaflets(leaflet_idx).extrusion_out = valve.extrusion_out; 
    end 
end

power_search_list = [2]; 
for leaflet_idx = 1:length(valve.leaflets)
    valve.leaflets(leaflet_idx).power_search_list = power_search_list; 
end 

% viscoelastic damping coefficients springs 
% eta, damping coeff here, is multiplied by the coefficient on the 
% associated spring 
% note that linear springs and collagen springs have vastly different constants 
% and these are tuned manually to make the dashpot constants equal order of magnitude
valve.eta_multiplier_linear   = 0; 
valve.eta_multiplier_collagen = 0; 

% valve_plot(valve); 
pause(.1); 

disp('Done with initialize.'); 


