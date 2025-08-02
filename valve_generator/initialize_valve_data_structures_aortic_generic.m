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

valve.rotate_identical_leaflets = true;

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

% does not place partition
valve.in_heart = true; 

fused_commissure = false;
cusp_fused = 0; 

% name 
if valve.in_heart
    valve.base_name = sprintf('aortic_no_partition_%d', N); 
    valve.extra_radius_hoops = 0.0; % adds points out the partition up to this amount 

    valve.tight_cylinder = true; 
    valve.z_extra_cylinder = 0.3; 
                                     
    % for normal_1
    % valve.initial_rotation_aortic = rotation_matrix_z(pi/4); 
    
    % for normal_3
    th = 2*pi/3; 
    valve.initial_translation_aortic = -0.05 * [cos(th); sin(th); 0]; 
    valve.initial_rotation_aortic = rotation_matrix_z(pi/3 + pi/12 + pi/48 + (cusp_fused * 2*pi/3));
    valve.transformation_vertex_file = 'aortic_annulus.vertex';
    
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
valve.num_copies = 3; 

% valve.copy_spring_weights = [1/2 1/4 1/4];

% respace on annulus in 3d 
% if false, spaced wrt theta 
valve.annulus_points_even_spacing = true; 

valve.use_annulus_flattened_pts = true; 

valve.annulus_to_comm = true; 

% from Khelil 2015 Surgical Anatomy of the Aortic Annulus Landmarks 
    % left half of nc leaflet plus reflections and normalization
valve.annulus_flattened_normalized = [ 
                   0   1.000000000000000
   0.002648301903414   0.829006318911995
   0.010063647338885   0.684870567517156
   0.034957485019150   0.522829189164894
   0.074152453295587   0.391226048708535
   0.113876981846794   0.298119648690557
   0.154661131477103   0.217546881096994
   0.224046541240633   0.126231145923352
   0.310381383503747   0.056401345909868
   0.382945155975022   0.020591113480958
   0.441207797850125   0.004476151367872
   0.499470439725229                   0
   0.500529560274771                   0
   0.558792202149874   0.004476151367872
   0.617054844024978   0.020591113480958
   0.689618616496253   0.056401345909868
   0.775953458759367   0.126231145923352
   0.845338868522897   0.217546881096994
   0.886123018153206   0.298119648690557
   0.925847546704413   0.391226048708535
   0.965042514980850   0.522829189164894
   0.989936352661115   0.684870567517156
   0.997351698096586   0.829006318911995
   1.000000000000000   1.000000000000000];

valve.normal_thicken = true; 
% nominal aortic valve thickness
valve.normal_thickness = 0.044 * (384/N); 

valve.extrusion_out = true;


valve.dirichlet_free_edge = false; 

valve.dirichlet_free_edge_with_ref_only = true; 



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

valve.p_final = 0 * MMHG_TO_CGS;  

% valve.dip_anterior_systole = true; 
% valve.r_dip = 0.75; 
% valve.total_angle_dip = pi; 

valve.L = 2.25; 

r_stj = 2.5/2; % 25 mm valve 
r_temp = 2.5/2; % vbr radius
hc = 0.5 * r_stj; 
h1 = 1.4 * r_stj - hc; 
r_commissure = r_stj; 
% place the post only if not using the full annulus geometry 
place_vertical_post = ~valve.annulus_to_comm;
valve.skeleton = get_skeleton_aortic_generic(r_temp, h1, hc, r_commissure, place_vertical_post); 
% valve.skeleton = get_skeleton_aortic_generic(); 
valve.r = valve.skeleton.r; 

valve.place_cylinder = true; 
valve.z_max_cylinder = (pi/3) * valve.r; 
valve.z_min_cylinder = 0.0; 
valve.n_layers_cylinder = 3; 


% Base constants, individual pieces are tuned relative to these values

% tension coefficients structure 

% pressure / tension coefficient ratio
% this tension coefficient is the maximum tension that a fiber can support
% valve.pressure_tension_ratio = 0.055; % 0.11 * 0.975; 
if valve.dirichlet_free_edge
    tension_coeffs.pressure_tension_ratio = 0.005; % 0.011; 
else 
    tension_coeffs.pressure_tension_ratio = 0.00477; 
end 

tension_coeffs.dec_tension_coeff_base = 20.0; 


% max tensions in leaflets 
tension_coeffs.alpha = 1.6;   % circumferential 
tension_coeffs.beta  = 0.055;   % radial

% decreasing tension coefficients 
tension_coeffs.c_circ_dec       = 4.8;  % circumferential 
tension_coeffs.c_rad_dec        = 1.53;  % radial

tension_coeffs.c_circ_dec_annulus = 1.8;        

% tension_coeffs.c_circ_dec_free_edge = 5.0;
tension_coeffs.c_circ_dec_free_edge_percentage = 0.0;

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
valve.eta_papillary_unscaled = 0.0; valve.target_papillary_unscaled/500; 

% if nonzero, linear springs of rest length with spacing between the layers 
% are placed with this value 
% final formula is multiplied by valve.tension_base  
valve.kappa_cross_layer_multipler = 10 * (384/N)^2 * 1e4 / 256^2;

% valve.k_bend_radial = [0 0 1e5 1e5] * 192/N;
valve.k_bend_radial = 1e4 * 192/N;
% valve.k_bend_radial_annulus = 1e2 * 192/N;
valve.k_bend_radial_free_edge = 0; 1e4 * 192/N;
valve.k_bend_radial_free_edge_percentage = 0; 
valve.k_bend_circ = 1e4 * 192/N; 
valve.k_bend_circ_free_edge = 0; 
valve.k_bend_circ_free_edge_percentage = 0;

valve.k_bend_cross_layer = 1e4 * 192/N;

if valve.in_heart
    
    dx = 5 /(N/4); 
    valve.ds = dx/2; %2*pi*valve.skeleton.r / N; 

    % if min radius lower than 2.5cm, increase ring thickness accordingly 
    thickness_cylinder = 0.3 + (2.5/2 - valve.r); 
    valve.n_layers_cylinder = ceil(thickness_cylinder/valve.ds) + 1; 

    h_scaffold_min = -0.05;

    h_top_scaffold_min = 0.05; 

    h_top_scaffold_max = valve.skeleton.normal_height; 

    h_top_scaffold_amplitude = h_top_scaffold_max - h_top_scaffold_min; 

    % value from least squares on pulm 
    % p = 3.095653361985474; 
    p = 100; 

    % function with unspecified power 
    % valve.z_max_cylinder = @(theta) h_top_scaffold_min * ones(size(theta))  +  h_top_scaffold_amplitude * abs(cos(theta)).^(p); 

    cos_power = @(theta) h_top_scaffold_min * ones(size(theta)) + (0.05 + h_top_scaffold_max) * abs(cos((3/2) * theta)).^(p); 

    annulus_min_fn = @(theta) h_top_scaffold_max * interp1(valve.annulus_flattened_normalized(:,1), valve.annulus_flattened_normalized(:,2), mod(theta,2*pi/3)/(2*pi/3), 'pchip'); 

    h_top_min_adjust = @(theta) h_top_scaffold_min * ones(size(theta));

    % extra_comm = @(theta) hc * abs(cos(theta)).^(p); 

    % valve.z_max_cylinder = @(theta) annulus_min_fn(theta) + extra_comm(theta) + h_top_min_adjust(theta);
    valve.z_max_cylinder = @(theta) max(cos_power(theta), annulus_min_fn(theta) + h_top_min_adjust(theta));

    % bottom flat 
    valve.z_min_cylinder = @(theta) h_scaffold_min * ones(size(theta)); 

    debug_plot = true; 
    if debug_plot
        figure; 
        th = linspace(0,2*pi,100000);
        plot(th,valve.z_min_cylinder(th))
        hold on 
        plot(th,valve.z_max_cylinder(th),'k')

        plot(th,cos_power(th))
        plot(th,annulus_min_fn(th))
%             plot(th,extra_comm(th))
%         f = @(theta)  0.28*ones(size(theta))  + (1.095 - 0.28)*0.5 * (cos(3*theta)+1);
%         plot(th, f(th)); 

        legend('bottom', 'top', 'cos_power', 'annulus_min')

%         plot(th, h_min_scaffold * ones(size(th))); 
%         plot(th, (h_min_scaffold + h_min_amplitude)*ones(size(th))); 
%         plot(th,  h_top_scaffold_min * ones(size(th))); 
%         plot(th,  h_top_scaffold_max * ones(size(th)));

        % heights from top of scaffold to minimum 
%        plot(angles, heights_theta,'k*')

        axis equal
    end 
%     dx = 5 /(N/4); 
%     valve.ds = dx/2; %2*pi*valve.skeleton.r / N; 
% 
%     thickness_cylinder = 0.15; 
%     valve.n_layers_cylinder = ceil(thickness_cylinder/valve.ds) + 1; 
%         
%     h_top_scaffold_min = 0.2; 
%     
%     h_top_scaffold_max = valve.skeleton.normal_height; 
%     
%     h_top_scaffold_amplitude = h_top_scaffold_max - h_top_scaffold_min; 
%     
%     % value from least squares on pulm 
%     % p = 3.095653361985474; 
%     p = 20; 
%     
%     % function with unspecified power 
%     valve.z_max_cylinder = @(theta) h_top_scaffold_min * ones(size(theta))  +  h_top_scaffold_amplitude * abs(cos((3/2)*theta)).^(p); 
%         
%     % bottom flat at zero 
%     valve.z_min_cylinder = @(theta) zeros(size(theta)); 
%     
%     debug_plot = false; 
%     if debug_plot
%         th = linspace(0,2*pi,1000);
%         plot(th,valve.z_min_cylinder(th))
%         hold on 
%         plot(th,valve.z_max_cylinder(th))
% 
% %         f = @(theta)  0.28*ones(size(theta))  + (1.095 - 0.28)*0.5 * (cos(3*theta)+1);
% %         plot(th, f(th)); 
%         
%         legend('bottom', 'top')
%         
% %         plot(th, h_min_scaffold * ones(size(th))); 
% %         plot(th, (h_min_scaffold + h_min_amplitude)*ones(size(th))); 
% %         plot(th,  h_top_scaffold_min * ones(size(th))); 
% %         plot(th,  h_top_scaffold_max * ones(size(th)));
%         
%         % heights from top of scaffold to minimum 
% %        plot(angles, heights_theta,'k*')
%         
%         axis equal
%         
%     end 
end 



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

[leaflet valve] = initialize_leaflet_aortic(name,                                ... 
                                            N,                                   ...
                                            tension_coeffs,                      ... 
                                            p_0,                                 ... 
                                            valve);  

valve.leaflets(1) = leaflet; 
valve.leaflets(1).extrusion_out = valve.extrusion_out;    

if fused_commissure
    valve.leaflets(1).fused_commissure = true; 
    valve.leaflets(1).fused_comm_idx    = 3; 
end 

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


