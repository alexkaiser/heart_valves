function [valve] = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential, decreasing_tension)
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
valve.max_it                = 4000; 
valve.max_it_continuation   = 2000; 

% Parameters for quick exit on line search 
valve.max_consecutive_fails = 5;  
valve.max_total_fails       = 40; 

if exist('attached', 'var') 
    valve.attached = attached; 
else 
    valve.attached = false; 
end 

% Valve skeleton parameters 
valve.r = 1.606587877768772; 

% original, supposedly diastolic but we have been using as systolic 
% valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683]; 
% valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348]; 

valve.left_papillary  = [ -0.972055648767080; -1.611924550017006; -2.990100960298683] + [0; 0; -0.0];
valve.right_papillary = [ -1.542417595752084;  1.611924550017006; -3.611254871967348] + [0; 0; -0.0]; 


% Places papillary attachments in linear interpolant between single point tips 

split_papillary = true; 

% vector pointing along line from left to right papillary 
l_to_r_papillary = (valve.right_papillary - valve.left_papillary); 
l_to_r_papillary = l_to_r_papillary / norm(l_to_r_papillary);


valve.papillary_radius = 0.25; 
 
valve.left_papillary_center  = valve.left_papillary  + valve.papillary_radius * l_to_r_papillary; 
valve.right_papillary_center = valve.right_papillary - valve.papillary_radius * l_to_r_papillary; 


diastolic_increment = [0; 0; 0.0]; 
valve.left_papillary_diastolic  = valve.left_papillary  + diastolic_increment; 
valve.right_papillary_diastolic = valve.right_papillary + diastolic_increment; 


valve.split_papillary = split_papillary; 
valve.radial_and_circumferential = true; 
valve.bead_slip = true; 
valve.leaflet_only = leaflet_only; 
valve.optimization = optimization; 
valve.repulsive_potential = repulsive_potential; 
valve.repulsive_power     = 1; 

if repulsive_potential
 
    % good total value (not including mesh parameters) at N=32
    repulsive_coeff_32 = 0.002238985441466; 
    
    repulsive_coeff_base = repulsive_coeff_32 * 32^2; 
    
    % scale so that when multiplied by above value gives the correct value 
    % valve.repulsive_coeff = repulsive_coeff_32 * 32^2; 
    
    valve.c_repulsive_circumferential = 2.0 * repulsive_coeff_base; 
    valve.c_repulsive_radial          = 6.0 * repulsive_coeff_base; 
    valve.c_repulsive_chordae         = 1.0 * repulsive_coeff_base; 
else 
    valve.repulsive_coeff  = 0.0; 
end 


valve.decreasing_tension = decreasing_tension; 

if decreasing_tension
 
    % good total value (not including mesh parameters) at N=32
    dec_tension_coeff_32 = 0.002238985441466; 
    
    dec_tension_coeff_base = dec_tension_coeff_32 * 32^2; 
    
    valve.c_dec_tension_circumferential = 1.0 * dec_tension_coeff_base; 
    valve.c_dec_tension_radial          = 1.0 * dec_tension_coeff_base; 
    valve.c_dec_tension_chordae         = 1.0 * dec_tension_coeff_base; 
else 
    valve.dec_tension  = 0.0; 
end 



% function pointers 
if attached 
    valve.diff_eqns = @difference_equations_bead_slip_attached; 
    valve.jacobian  = @build_jacobian_bead_slip_attached; 
    
    if leaflet_only
        error('leaflet_only not implemented for attached')
    end 
        
else 
    if leaflet_only
        valve.diff_eqns = @difference_equations_bead_slip_leaflet_only; 
        valve.jacobian  = @build_jacobian_bead_slip_leaflet_only;
    else 
        valve.diff_eqns = @difference_equations_bead_slip; 
        valve.jacobian  = @build_jacobian_bead_slip;
        
    end 
end 

% general solve parameters

% name 
valve.base_name = sprintf('mitral_tree_%d', N); 

% box width 
valve.L = 2.5; 

% pressure / tension coefficient ratio
% this tension coefficient is the maximum tension that a fiber can support
valve.pressure_tension_ratio = 0.15; % 0.11 * 0.975; 

% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
valve.refinement = N/32.0; 

MMHG_TO_CGS = 1333.22368;
valve.p_physical = 110 * MMHG_TO_CGS; 

% scaling for target points 
valve.target_multiplier = 40/128; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
valve.n_lagrangian_tracers = 0; 

% Uses configuration of X 
valve.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
valve.num_copies = 1; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
valve.collagen_constitutive = true; 


% anterior leaflet data structure 
reflect_x = false; 
total_angle_anterior = 5*pi/6; 
angles_anterior = [-total_angle_anterior/2, total_angle_anterior/2]; 

% Radial and circumferential fibers 
% Or diagonally oriented fibers 
radial_and_circumferential = true; 

% base constant for tensions 
valve.tension_base = valve.p_physical / valve.pressure_tension_ratio; 

% physical units create a scalar multiple of the old 
% this multiple is large number, so we want to scale the old tolerance accordingly 
% 8.3326e-04 is a good number here
valve.tol_global = 1e-3; 

% pressure on each leaflet is constant, and always a physical pressure 
p_0 = -valve.p_physical; 


% Spring constants in two directions 
tension_base_anterior = valve.tension_base; 
alpha_anterior    = 1.0 * tension_base_anterior;  % circumferential 
beta_anterior     = 1.0 * tension_base_anterior;  % radial


% Places this many extra fibers from ring to ring 
% Must be 0 <= N_ring_to_ring <= (N/2)
ring_to_ring_anterior_range = 0; %[(N/8), (3*N/8)];


% Add energy function for zero pressure case 
if (p_0 == 0.0) && (~leaflet_only)
    valve.energy = @energy_bead_slip; 
end 


% Chordae parameters 

n_trees_anterior = 4; 

% total force in leaves of each tree 
% good basic value = 1.8 * 0.5 * (alpha + beta) = 1.8 * tension_base
k_0_1_anterior = 0.8 * 2.0 * tension_base_anterior / n_trees_anterior; 

% try turning up outside values 
k_0_1_anterior = k_0_1_anterior * [1 1 1 1]; 

% constant tension at the root of the tree 
% this is determined by hand tuning k_multiplier at coarse resolution 
% then taking the k_root 
% base good value, old program units  
% k_root = 1.889568000000001e+01 / 32; 
% adjust accordingly
% note that 1.889568000000001e+01 / 32 = 0.5905
k_root_anterior = 0.9 * 2.0 * (1.889568000000001e+01 / 32) * tension_base_anterior / n_trees_anterior; 

k_root_anterior = k_root_anterior * [1 1 1 1]; 

% controls initial guess tree vertex placement 
tree_frac = 0.5;


papillary_anterior = zeros(3,n_trees_anterior); 

n_points = n_trees_anterior/2; 

right_papillary_range = 1:(n_trees_anterior/2); 
left_papillary_range  = right_papillary_range + (n_trees_anterior/2);

papillary_anterior(:,right_papillary_range) = get_papillary_coords(valve.left_papillary_center,  valve.papillary_radius, n_points,  0*pi/4,    pi/4); 
papillary_anterior(:,left_papillary_range)  = get_papillary_coords(valve.right_papillary_center, valve.papillary_radius, n_points,   -pi/4, -0*pi/4); 

n_leaves_anterior  = N/n_trees_anterior * ones(n_trees_anterior, 1); 
tree_direction_anterior = [-1; -1; 1; 1];

left_papillary_anterior_diastolic  = valve.left_papillary_diastolic; 
right_papillary_anterior_diastolic = valve.right_papillary_diastolic;


valve.anterior = initialize_leaflet_bead_slip(N,                        ... 
                                    reflect_x,                          ... 
                                    angles_anterior,                    ...    
                                    valve.r,                            ... 
                                    papillary_anterior,                 ... 
                                    n_leaves_anterior,                  ...
                                    tree_direction_anterior,            ...
                                    left_papillary_anterior_diastolic,  ... 
                                    right_papillary_anterior_diastolic, ... 
                                    radial_and_circumferential,         ...  
                                    alpha_anterior,                     ... 
                                    beta_anterior,                      ...
                                    p_0,                                ... 
                                    k_0_1_anterior,                     ... 
                                    k_root_anterior,                    ... 
                                    tree_frac,                          ... 
                                    leaflet_only,                       ...
                                    ring_to_ring_anterior_range,        ...
                                    valve);  

                   


    
total_angle_posterior = 7*pi/6; 
tension_base_posterior = valve.tension_base; 
k_0_1_posterior  = 0.2 * tension_base_posterior; 
k_root_posterior = 0.8 * 2.0 * (1.889568000000001e+01 / 32) * tension_base_posterior; 
ring_to_ring_posterior_range = 0; %[(3*N/16), (3*N/8)];
     
    
angles_posterior = [pi - total_angle_posterior/2, pi + total_angle_posterior/2]; 
reflect_x = false;  

alpha_posterior    = 1.0 * tension_base_posterior;  % circumferential 
beta_posterior     = 1.0 * tension_base_posterior;  % radial


n_trees_posterior = 8; 

k_root_posterior = k_root_posterior / n_trees_posterior; 

N_tree = N/n_trees_posterior; 



papillary_posterior = zeros(3,n_trees_posterior); 

n_points = n_trees_posterior/2; 

right_papillary_range = 1:(n_trees_posterior/2); 
left_papillary_range  = right_papillary_range + (n_trees_posterior/2); 

papillary_posterior(:,right_papillary_range) = get_papillary_coords(valve.right_papillary_center, valve.papillary_radius, n_points,    pi/4,  5*pi/4); 
papillary_posterior(:,left_papillary_range)  = get_papillary_coords(valve.left_papillary_center,  valve.papillary_radius, n_points, -5*pi/4,   -pi/4);

% this arrangement is very touchy, doesn't converge  
% N * [-1, 1, -1, 1, -1, 1] .* [1/8, 1/8, 1/4, 1/4, 1/8, 1/8]; 

% this is generally pretty good 
n_leaves_posterior = N/n_trees_posterior * ones(n_trees_posterior, 1); 
tree_direction_posterior = [-1; -1; -1; -1; 1; 1; 1; 1]; 


% n_leaves_and_direction_posterior = N * [1/8, 1/8, 1/4, 1/4, 1/8, 1/8] .* [-1, -1, -1, 1, 1, 1]; 

left_papillary_posterior_diastolic  = valve.left_papillary_diastolic; 
right_papillary_posterior_diastolic = valve.right_papillary_diastolic;


valve.posterior = initialize_leaflet_bead_slip(N,                    ...
                                reflect_x,                           ... 
                                angles_posterior,                    ...    
                                valve.r,                             ... 
                                papillary_posterior,                 ... 
                                n_leaves_posterior,                  ...
                                tree_direction_posterior,            ...
                                left_papillary_posterior_diastolic,  ...
                                right_papillary_posterior_diastolic, ...
                                radial_and_circumferential,          ...  
                                alpha_posterior,                     ... 
                                beta_posterior,                      ... 
                                p_0,                                 ... 
                                k_0_1_posterior,                     ... 
                                k_root_posterior,                    ... 
                                tree_frac,                           ... 
                                leaflet_only,                        ...
                                ring_to_ring_posterior_range,        ...
                                valve);  

    
    

valve_plot(valve); 

disp('Done with initialize.'); 


