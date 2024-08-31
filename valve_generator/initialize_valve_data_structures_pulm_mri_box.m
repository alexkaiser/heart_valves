function [valve] = initialize_valve_data_structures_pulm_mri_box(N)
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

mri_box = true;

% name 
if valve.in_heart
    valve.base_name = sprintf('aortic_no_partition_%d', N); 
    valve.extra_radius_hoops = 0.0; % adds points out the partition up to this amount 

    valve.tight_cylinder = true; 
    valve.z_extra_cylinder = 0.28; 
                                     
    if mri_box
        valve.initial_translation_aortic = [0; 0; -0.2]; % 2mm down from origin, translation applied first 
        valve.initial_rotation_aortic = rotation_matrix_y(-pi/2) * rotation_matrix_z(pi/6); 
    else 
        valve.initial_rotation_aortic = rotation_matrix_z(pi/4); 
        valve.transformation_vertex_file = 'aortic_annulus.vertex';
    end 
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

valve.normal_thicken = true; 
% nominal aortic valve thickness
valve.normal_thickness = 0.044; 

valve.center_extrusion = true; 

valve.extrusion_out = false; 

valve.pre_extrude = false; 

% add flags to spring files 
% to view and output with a stride 

valve.output.leaflets       = [1;1;1]; 
valve.output.stride_leaflet = max(1,N/128); 
valve.output.mesh           = [1;0;0]; 
valve.output.cartesian_mesh = [0;0;0]; 
valve.output.stride_mesh    = N/32; 


valve.dirichlet_free_edge = false; 

valve.dirichlet_free_edge_with_ref_only = true; 

pinch_commissure = true; 
N_points_half_free_edge = (N/3)/2; 
N_to_pinch = N_points_half_free_edge/4;  

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
valve.p_physical = 30 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

valve.p_final = 0 * MMHG_TO_CGS;  

% valve.dip_anterior_systole = true; 
% valve.r_dip = 0.75; 
% valve.total_angle_dip = pi; 

valve.L = 2.25; 


r = 1.0; 
normal_height = 1.1; % 0.845; 
hc = 0.3; % 1 mm of commissure attachment (nearly zero)
h1 = normal_height - hc; 
valve.skeleton = get_skeleton_aortic_generic(r, h1, hc);

valve.r = valve.skeleton.r; 

% little nub at top of valve 
r_subtract_nub = 0.15; 
valve.skeleton.r_of_z = @(z) r .* ones(size(z)) - (abs(z - 1.2) < .1) .* r_subtract_nub .* cos( (pi/2)*(z - 1.2)/.1 ); 
                         
                         
r_of_z_debug = true; 
if r_of_z_debug 
    z_range = linspace(-.2, 1.3, 1000); 
    figure; 
    plot(valve.skeleton.r_of_z(z_range), z_range); 
end

                         

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
tension_coeffs.c_circ_dec       = 5.65;  % circumferential 
tension_coeffs.c_rad_dec        = 7.0;  % radial

tension_coeffs.c_circ_dec_annulus = 2.1;        

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
valve.kappa_cross_layer_multipler = (384/N)^2 * (1e4 / 256^2); 

% valve.k_bend_radial = [0 0 1e5 1e5] * 192/N;
valve.k_bend_radial = 0; 1e2 * 192/N;
% valve.k_bend_radial_annulus = 1e2 * 192/N;
valve.k_bend_radial_free_edge = 0; 1e4 * 192/N;
valve.k_bend_radial_free_edge_percentage = 0; 
valve.k_bend_circ = 0; 
valve.k_bend_circ_free_edge = 0; 
valve.k_bend_circ_free_edge_percentage = 0;

valve.k_bend_cross_layer = 0;

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
dx = 5 /(N/4); 
valve.ds = dx/2; %2*pi*valve.skeleton.r / N; 

if mri_box
    thickness_cylinder = 0.3 + r_subtract_nub; 
    valve.n_layers_cylinder = ceil(thickness_cylinder/valve.ds) + 1; 
    
    valve.r_max_cylinder = 1.3; 
    
    % top scaffold heights relative to table 
    h_top_scaffold_min_table = 0.65; 
    h_top_scaffold_max_table  = 1.5; 
    
    % bottom of scaffold goes up by this much at the supports 
    h_min_amplitude = 0.2; 
    h_min_scaffold_table = 0; % zero by definition, since this rests on the table  
    
    % origin is placed at this height relative to table 
    % this is not required to be equal
    origin_height_table = 0.2; 
    
    h_min_scaffold     = h_min_scaffold_table     - origin_height_table; 
    h_top_scaffold_min = h_top_scaffold_min_table - origin_height_table; 
    h_top_scaffold_max = h_top_scaffold_max_table - origin_height_table; 
    
    h_top_scaffold_amplitude = h_top_scaffold_max - h_top_scaffold_min; 
    
    % top gets cosine to a power 
    % least squares fit this 
    % valve.z_max_cylinder = @(theta) h_top_scaffold_min * ones(size(theta))  +  h_top_scaffold_amplitude * abs(cos((3/2)*theta)).^4; 
    
    % function with unspecified power 
    z_max_tmp = @(p,theta) h_top_scaffold_min * ones(size(theta))  +  h_top_scaffold_amplitude * abs(cos((3/2)*theta)).^(p); 
    
    % heights measured 
    heights_theta = [1.5 1.15 0.8 0.65] - origin_height_table;
    angles = [0:3] * 2*pi/18;
    
    % least squares fit it 
    myfittype = fittype(z_max_tmp,...
    'dependent',{'h_tmp'},'independent',{'theta'},...
    'coefficients',{'p'}); 
    
    myfit = fit(angles',heights_theta',myfittype); 
    
    % evaluate at least squares value of p 
    valve.z_max_cylinder = @(theta) z_max_tmp(myfit.p, theta); 
    
    % bottom just looks like a cosine 
    valve.z_min_cylinder = @(theta) h_min_scaffold * ones(size(theta)) +  h_min_amplitude * 0.5*(cos(3*theta)+1); 
    
    debug_plot = true; 
    if debug_plot
        fig = figure; 
        th = linspace(0,2*pi,1000);
        plot(th,valve.z_min_cylinder(th))
        hold on 
        plot(th,valve.z_max_cylinder(th))

%         f = @(theta)  0.28*ones(size(theta))  + (1.095 - 0.28)*0.5 * (cos(3*theta)+1);
%         plot(th, f(th)); 
        
        legend('bottom', 'top')
        
        plot(th, h_min_scaffold * ones(size(th))); 
        plot(th, (h_min_scaffold + h_min_amplitude)*ones(size(th))); 
        plot(th,  h_top_scaffold_min * ones(size(th))); 
        plot(th,  h_top_scaffold_max * ones(size(th)));
        
        % heights from top of scaffold to minimum 
        plot(angles, heights_theta,'k*')

        axis equal
        printfig(fig, "annulus morphology")
    end 
    
    
    
    
    
end 

[leaflet valve] = initialize_leaflet_aortic(name,                                ... 
                                            N,                                   ...
                                            tension_coeffs,                      ... 
                                            p_0,                                 ... 
                                            valve);  

valve.leaflets(1) = leaflet; 
    
if pinch_commissure
    valve.leaflets(1).pinch_commissure = true; 
    valve.leaflets(1).N_to_pinch       = N_to_pinch; 
end 

if isfield(valve, 'extrusion_out') 
    valve.leaflets(1).extrusion_out = valve.extrusion_out; 
end

if isfield(valve, 'pre_extrude')
    valve.leaflets(1).pre_extrude = valve.pre_extrude; 
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


