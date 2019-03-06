% check_jacobian_taylor_series
% 
% checks jacobian is actually an approximation to the derivative
% to the expected order using a Taylor series 
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


% reset stream for consistent results 

N = 8; 

% Initialize structures  
attached = false; 
leaflet_only = false; 
valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only); 

rand('twister',76599)

epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 

leaflet = valve.anterior; 

% eval the difference eqns on the perturbation 
[F_anterior F_chordae_left F_chordae_right] = difference_equations_bead_slip(leaflet); 

F_linearized = linearize_internal_points(leaflet, F_anterior, F_chordae_left, F_chordae_right);

% jacobian does not change 
J = build_jacobian_bead_slip(leaflet); 

fig = figure; 
spy(J, 'k'); 
title('Jacobian nonzero structure')

j_max       = leaflet.j_max; 
k_max       = leaflet.k_max; 
is_internal = leaflet.is_internal; 

% perturbation also does not change 
Z = zeros(size(leaflet.X)); 
for j=1:j_max
    for k=1:k_max
        if is_internal(j,k)
            Z(:,j,k) = rand(3,1);  
        end 
    end 
end 

leaflet_Z   = leaflet; 
leaflet_Z.X = Z; 
leaflet_Z.chordae.C_left  = rand(size(leaflet_Z.chordae.C_left)); 
leaflet_Z.chordae.C_right = rand(size(leaflet_Z.chordae.C_right)); 
Z_linearized = linearize_internal_points(leaflet_Z, leaflet_Z.X, leaflet_Z.chordae.C_left, leaflet_Z.chordae.C_right); 



fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    leaflet_perturbation   = leaflet;  
    leaflet_perturbation.X = leaflet.X + ep*leaflet_Z.X; 
    
    leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
    leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

    % eval the difference eqns on the perturbation 
    [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip(leaflet_perturbation); 
    F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

    
    errors(i) = norm(F_perturbed_linearized - F_linearized - ep*J*Z_linearized, 2); 
    
    fprintf('%e\t | %e \n', ep, errors(i)); 

end 

fprintf('\n\n\n\n'); 

figure; 
loglog(epsilon_vals, errors, '-*'); 
hold on 
loglog(epsilon_vals, epsilon_vals.^2, '--'); 

legend('error', 'eps^2')


figure; 

% leaflet part 
for k=1:k_max
    for j=1:j_max
        if is_internal(j,k)

            j
            k

            errors = zeros(size(epsilon_vals)); 

            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                leaflet_perturbation   = leaflet;  
                leaflet_perturbation.X = leaflet.X + ep * leaflet_Z.X; 

                leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
                leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

                % eval the difference eqns on the perturbation 
                [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip(leaflet_perturbation); 
                F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = leaflet_perturbation.linear_idx_offset(j,k) + (1:3);                     
                errors(i) = norm(diffs(range)); 

                fprintf('%e\t | %e \n', ep, errors(i)); 

            end 

            fprintf('\n\n\n\n'); 

            loglog(epsilon_vals, errors, '-*'); 
            hold on 
            loglog(epsilon_vals, epsilon_vals.^2, '--'); 

            legend('error', 'eps^2')

        end 
    end 
end 


% chordae part if included 
[m N_chordae] = size(leaflet.chordae.C_left); 
total_internal = 3*sum(is_internal(:)); 

for left_side = [true false]
    for i=1:N_chordae

        left_side
        i

        errors = zeros(size(epsilon_vals)); 


        for ep_idx = 1:length(epsilon_vals)

            ep = epsilon_vals(ep_idx); 

            % make a new structure for the perturbation 
            leaflet_perturbation   = leaflet;  
            leaflet_perturbation.X = leaflet.X + ep * leaflet_Z.X; 

            leaflet_perturbation.chordae.C_left  = leaflet.chordae.C_left  + ep*leaflet_Z.chordae.C_left; 
            leaflet_perturbation.chordae.C_right = leaflet.chordae.C_right + ep*leaflet_Z.chordae.C_right; 

            % eval the difference eqns on the perturbation 
            [F_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip(leaflet_perturbation); 
            F_perturbed_linearized = linearize_internal_points(leaflet_perturbation, F_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

            diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

            range = range_chordae(total_internal, N_chordae, i, left_side); 
            errors(ep_idx) = norm(diffs(range)); 

            fprintf('%e\t | %e \n', ep, errors(ep_idx)); 

        end 

        fprintf('\n\n\n\n'); 

        loglog(epsilon_vals, errors, '-*'); 
        hold on 
        loglog(epsilon_vals, epsilon_vals.^2, '--'); 

        legend('error', 'eps^2')


    end 
end 
 












