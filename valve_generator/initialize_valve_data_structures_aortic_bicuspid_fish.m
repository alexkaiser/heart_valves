function [valve] = initialize_valve_data_structures_aortic_bicuspid_fish(N)
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

valve.targets_for_bcs = false; 
valve.targets_for_bcs_ref_only = false; 

% general solve parameters
name = 'aortic'; 
valve.name = name; 

variety= 'bicuspid'; 
valve.variety = variety; 

% does not place partition
valve.in_heart = true; 
graft_tester_geometry = true; 


% name 
if valve.in_heart
    valve.base_name = sprintf('aortic_fish_%d', N); 
    valve.extra_radius_hoops = 0.0; % adds points out the partition up to this amount 

    if ~graft_tester_geometry    
        valve.tight_cylinder = true; 
        valve.z_extra_cylinder = 0.3; 
    end 
else 
    valve.base_name = sprintf('aortic_fish_box_%d', N); 
end 
MMHG_TO_CGS     = 1333.22368;


% % number of lagrangian tracers in each dimension 
% % arranged in a mesh near the origin
% % z direction is doubled 
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


% add flags to spring files 
% to view and output with a stride 

valve.output.leaflets       = [1;1;1]; 
valve.output.stride_leaflet = max(1,N/128); 
valve.output.mesh           = [1;0;0]; 
valve.output.cartesian_mesh = [0;0;0]; 
valve.output.stride_mesh    = N/32; 


valve.dirichlet_free_edge = false; 

valve.dirichlet_free_edge_with_ref_only = true; 


% provides a bending resistance for the final solve 
% for initial conditions 
% this is just to get a reasonable initial condition 
valve.k_bend_radial_ref_only = 0; 1e-14 * (96/N); 

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

% mean fwd pressure diff 
valve.p_physical = 0.9250 * MMHG_TO_CGS; 

% Pressure on each leaflet is constant, negative since normal is outward facing 
p_0 = -valve.p_physical; 

valve.p_final = 0.0 * MMHG_TO_CGS;  


valve.dirichlet_free_edge_comm_ref_only = false; 
valve.n_fixed_comm = max(1, floor(8*N/192));  
valve.p_final_fixed_comm = 0.1 * MMHG_TO_CGS;  


% valve.dip_anterior_systole = true; 
% valve.r_dip = 0.75; 
% valve.total_angle_dip = pi; 

% scale human box length down by 100 
valve.L = 2.25 / 100; 

% scan points 
% valve.skeleton = get_skeleton_fish(); 

% 300 microns 
r_temp = 363e-4 / 2; 
hc = 0.1 * r_temp; 
h1 = 285e-4 - hc; 
valve.skeleton = get_skeleton_aortic_generic(r_temp, h1, hc); 
valve.r = valve.skeleton.r; 
% 
valve.place_cylinder = true; 
valve.z_max_cylinder = 1.4 * valve.r; 
valve.z_min_cylinder = 0.0; 


% proportional takedown from 25 mm diameter valve 
% distance_scaling = 1e-2; 
distance_scaling = (valve.r * 2) / 2.5; 

valve.normal_thicken = true; 
% nominal aortic valve thickness
valve.normal_thickness = 0.044 * distance_scaling; 

valve.extrusion_out = true; 

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
tension_coeffs.c_circ_dec       = 5;  % circumferential 
tension_coeffs.c_rad_dec        = 2;  % radial

tension_coeffs.c_circ_dec_annulus = 1.31;

% tension_coeffs.c_circ_dec_free_edge = 5.0;
tension_coeffs.c_circ_dec_free_edge_percentage = 0.0;

% scaling for target points 
% note that this does not include copies 
% and scaling for copies is handled by the output routine 

% scales for by mesh width for consistant total mesh force on ring 
valve.target_net_unscaled       = 1e-3 * (8 / valve.N) * (192/N); 

% does not scale since total number of points is constant 
valve.target_papillary_unscaled = 0; 

% viscoelastic damping coefficients for net, does not include copies 
valve.eta_net_unscaled = 0; % 1e-5 * valve.target_net_unscaled; 

% viscoelastic damping coefficients for root attachments, does not include copies  
valve.eta_papillary_unscaled = 0.0; valve.target_papillary_unscaled/500; 

% if nonzero, linear springs of rest length with spacing between the layers 
% are placed with this value 
% final formula is multiplied by valve.tension_base  
valve.kappa_cross_layer_multipler = 1e-3 * 2 * (384/N)^2 * 1e4 / 256^2;


if valve.in_heart 

    if graft_tester_geometry        
    
        distance_scaling_fluid = 1e-2; 
        
        dx = distance_scaling_fluid * 0.1 * (192/N);
        NZ = 128 * (N/192); 
        valve.z_min_cylinder = -3 * distance_scaling_fluid; 
        valve.z_max_cylinder = valve.z_min_cylinder + dx * (NZ - 1); 
        
        % update r_of_z for extender 
        
        extender_extra_rad = 0.5 * distance_scaling_fluid; 
        extender_length = 1.5 * distance_scaling_fluid; 
        
        valve.skeleton.r_of_z = @(z) valve.skeleton.r .* ones(size(z)) + ...
                                     (z <= (valve.z_min_cylinder + extender_length)) .* ... % mask for bottom portion 
                                     extender_extra_rad .* 0.5 .* (cos(pi * (z - valve.z_min_cylinder)/(extender_length)) + 1.0); 

        
        debug_extender_plot = false; 
        if debug_extender_plot 
            z = valve.z_min_cylinder:0.0001:valve.z_max_cylinder;             
            figure; 
            plot(valve.skeleton.r_of_z(z), z)
            ylabel('z')
            xlabel('r_of_z')
            xlim([0 3])
            axis equal 
            
        end
    else 
        error('not implemented');         
    end
end 





% Approximate Lagrangian mesh spacing at ring 
% Used for later splitting of springs 
% If any spring is placed at more than double this length an extra vertex is placed
valve.ds = 2*pi*valve.skeleton.r / N; 

valve.k_rel_unscaled = valve.ds * 1e-3; 

[leaflet valve] = initialize_leaflet_aortic(name,                                ... 
                                            N,                                   ...
                                            tension_coeffs,                      ... 
                                            p_0,                                 ... 
                                            valve);  

valve.leaflets(1) = leaflet; 
    
if isfield(valve, 'extrusion_out') 
    valve.leaflets(1).extrusion_out = valve.extrusion_out; 
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


