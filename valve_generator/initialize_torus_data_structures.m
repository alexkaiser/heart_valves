function [torus] = initialize_torus_data_structures(N, repulsive_potential)
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
torus.N = N; 
torus.j_max = N; 
torus.k_max = N; 
torus.tol_global = 1e-10;
torus.max_it = 4000; 

torus.bead_slip = true; 
torus.leaflet_only = true; 

torus.repulsive_potential = false; 

if repulsive_potential
    error('repulsive potential not implemented for torus')
end 

torus.repulsive_power     = 1; 

% general mesh parameters 
torus.du = 1/N; 
torus.dv = 1/N; 

% fiber spacing util parameter 
dt = 1/N; 

% coefficient has units of 1/L, then gets squared 
c_0 = N;  
torus.repulsive_coeff  = 1e-7 * c_0^2; 

torus.diff_eqns = @difference_equations_torus; 
torus.jacobian  = @build_jacobian_torus; 
    
% major and minor radius of torus 
torus.r = 1; 
torus.R = 2;

% number of wraps in major and minor directions 

% minor direction may wrap zero times 
torus.n_wraps = 1; 

% major direction must wrap at least once 
torus.m_wraps = 2; 


% general solve parameters

% name 
torus.base_name = sprintf('torus_%d', N); 

% box width 
torus.L = 2.5; 

% pressure / spring constant ratio  
% ratio 6 is for N=32
% ratio = 6 seems to make everything very stiff 
% turn down by order of magnitude, see if it helps 
torus.pressure_tension_ratio = 1.5; 


% original spring constants were for N = 32 debug width
% spring constants get multiplied by 32/N, so they are halfed if N==64
% use this refintement number accordingly 
torus.refinement = N/32.0; 

torus.p_physical = 100; 

% scaling for target points 
torus.target_multiplier = 40; 

% number of lagrangian tracers in each dimension 
% arranged in a mesh near the origin
% z direction is doubled 
torus.n_lagrangian_tracers = 8; 

% Uses configuration of X 
torus.X_config_is_reference = true; 

% places this many exact copies of the leaflet downward in z 
% spring constants are all reduced by num_copies 
% spacing is always half a mesh width 
torus.num_copies = 3; 

% Uses collagen spring function implemented in IBAMR 
% Spring constants are different here 
torus.collagen_springs_leaflet = false; 


% Spring constants in two directions 
torus.alpha    =  1.0;   
torus.beta     =  1.0;  

% parametrization gives outward normal 
% positive p_0 gives 
torus.p_0      =  0.1;  



tor = @(u,v) [ cos(v) .* (torus.R + torus.r * cos(u)); ... 
               sin(v) .* (torus.R + torus.r * cos(u)); ...  
                                    torus.r * sin(u)]; 


torus.preimage = zeros(2,N,N);                                 
                                
torus.X = zeros(3,N,N);


for k=1:N
    
    % initial u conditions go from zero 
    u_0 = 2*pi * (k-1) * dt; 
    
    for j=1:N
        
        u = 2*pi*torus.n_wraps * dt * (j-1) + u_0; 
        v = 2*pi*torus.m_wraps * dt * (j-1); 
                
        torus.preimage(:,j,k) = [u; v];  
        
        torus.X(:,j,k) = tor(u,v); 
    end 
end 

fig = figure; 
fig = torus_plot(torus, fig); 


'done with initialize'




fig = figure; 
hold on 

for k=1:N
    plot(torus.preimage(1,:,k), torus.preimage(2,:,k), 'o'); 
end 
axis equal 
title('preimage')


% set util arrays 
torus.is_internal       =  ones(torus.j_max, torus.k_max); 
torus.is_bc             = zeros(torus.j_max, torus.k_max); 
torus.linear_idx_offset = zeros(torus.j_max, torus.k_max); 
torus.point_idx_with_bc = zeros(torus.j_max, torus.k_max); 

count = 0; 
for k=1:torus.k_max
    for j=1:torus.j_max
        if torus.is_internal(j,k)
            torus.linear_idx_offset(j,k) = count; 
            count = count + 3; 
        end 
    end 
end

count = 0;
for k=1:torus.k_max
    for j=1:torus.j_max
        if torus.is_internal(j,k) || torus.is_bc(j,k)
            torus.point_idx_with_bc(j,k) = count; 
            count = count + 1; 
        end 
    end 
end

torus.chordae.C_left = []; 
torus.chordae.C_right = []; 


data_movement_check = false; 
if data_movement_check

    'version with no movement, including b.c.s'
    torus.X 

    % 'linear order on internal points'
    X_linearized = linearize_internal_points(torus, torus.X) 

    
    'moved back to 3d ordering'
    torus = internal_points_to_2d(X_linearized, torus); 
    torus.X 

end 








