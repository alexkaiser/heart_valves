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
radial       = true; 
bead_slip    = true; 
attached     = true; 
leaflet_only = false; 
optimization = false; 
repulsive_potential = false;

valve = initialize_valve_data_structures_radial_bead_slip(N, attached, leaflet_only, optimization, repulsive_potential); 


rand('twister',76599)

epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 



% eval the difference eqns on the perturbation 
[F_anterior F_posterior F_chordae_left F_chordae_right] = difference_equations_bead_slip_attached(valve); 
F_linearized = linearize_internal_points_bead_slip_attached(valve, F_anterior, F_posterior, F_chordae_left, F_chordae_right);

% jacobian does not change 
J = build_jacobian_bead_slip_attached(valve); 

figure; 
spy(J); 
title('jacobian nonzero structure in jacobian tester')

anterior = valve.anterior; 
posterior = valve.posterior; 


j_max                = anterior.j_max; 
k_max                = anterior.k_max; 
is_internal_anterior = anterior.is_internal; 

% perturbation also does not change 
Z_anterior = zeros(size(anterior.X)); 
for j=1:j_max
    for k=1:k_max
        if is_internal_anterior(j,k)
            Z_anterior(:,j,k) = rand(3,1);  
        end 
    end 
end 


j_max                = posterior.j_max; 
k_max                = posterior.k_max; 
is_internal_posterior = posterior.is_internal; 

% perturbation also does not change 
Z_posterior = zeros(size(posterior.X)); 
for j=1:j_max
    for k=1:k_max
        if is_internal_posterior(j,k)
            Z_posterior(:,j,k) = rand(3,1);  
        end 
    end 
end 


valve_Z   = valve; 
valve_Z.anterior.X = Z_anterior; 
valve_Z.posterior.X = Z_posterior; 
valve_Z.anterior.chordae.C_left  = rand(size(valve_Z.anterior.chordae.C_left)); 
valve_Z.anterior.chordae.C_right = rand(size(valve_Z.anterior.chordae.C_right));

Z_linearized = linearize_internal_points_bead_slip_attached(valve_Z, valve_Z.anterior.X, valve_Z.posterior.X, valve_Z.anterior.chordae.C_left, valve_Z.anterior.chordae.C_right); 


fprintf('eps\t | taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    valve_perturbation   = valve;  
    valve_perturbation.anterior.X  = valve.anterior.X  + ep*valve_Z.anterior.X; 
    valve_perturbation.posterior.X = valve.posterior.X + ep*valve_Z.posterior.X; 
    
    valve_perturbation.anterior.chordae.C_left  = valve.anterior.chordae.C_left  + ep*valve_Z.anterior.chordae.C_left; 
    valve_perturbation.anterior.chordae.C_right = valve.anterior.chordae.C_right + ep*valve_Z.anterior.chordae.C_right; 

    % eval the difference eqns on the perturbation 
    [F_anterior_perturbed F_posterior_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip_attached(valve_perturbation); 
        
    F_perturbed_linearized = linearize_internal_points_bead_slip_attached(valve, F_anterior_perturbed, F_posterior_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 
    
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
        if is_internal_anterior(j,k)

            'anterior'
            j
            k
            

            errors = zeros(size(epsilon_vals)); 

            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                valve_perturbation   = valve;  
                valve_perturbation.anterior.X  = valve.anterior.X  + ep*valve_Z.anterior.X; 
                valve_perturbation.posterior.X = valve.posterior.X + ep*valve_Z.posterior.X; 

                valve_perturbation.anterior.chordae.C_left  = valve.anterior.chordae.C_left  + ep*valve_Z.anterior.chordae.C_left; 
                valve_perturbation.anterior.chordae.C_right = valve.anterior.chordae.C_right + ep*valve_Z.anterior.chordae.C_right; 

                % eval the difference eqns on the perturbation 
                [F_anterior_perturbed F_posterior_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip_attached(valve_perturbation); 

                F_perturbed_linearized = linearize_internal_points_bead_slip_attached(valve, F_anterior_perturbed, F_posterior_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = valve.anterior.linear_idx_offset(j,k) + (1:3);                     
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



% leaflet part 
for k=1:k_max
    for j=1:j_max
        if is_internal_posterior(j,k)

            'posterior'
            j
            k
            

            errors = zeros(size(epsilon_vals)); 

            for i = 1:length(epsilon_vals)

                ep = epsilon_vals(i); 

                % make a new structure for the perturbation 
                valve_perturbation   = valve;  
                valve_perturbation.anterior.X  = valve.anterior.X  + ep*valve_Z.anterior.X; 
                valve_perturbation.posterior.X = valve.posterior.X + ep*valve_Z.posterior.X; 

                valve_perturbation.anterior.chordae.C_left  = valve.anterior.chordae.C_left  + ep*valve_Z.anterior.chordae.C_left; 
                valve_perturbation.anterior.chordae.C_right = valve.anterior.chordae.C_right + ep*valve_Z.anterior.chordae.C_right; 

                % eval the difference eqns on the perturbation 
                [F_anterior_perturbed F_posterior_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip_attached(valve_perturbation); 

                F_perturbed_linearized = linearize_internal_points_bead_slip_attached(valve, F_anterior_perturbed, F_posterior_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

                diffs = F_perturbed_linearized - F_linearized - ep*J*Z_linearized; 

                range = valve.posterior.linear_idx_offset(j,k) + (1:3);                     
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
[m N_chordae] = size(valve.anterior.chordae.C_left); 
total_internal = 3*(sum(is_internal_anterior(:)) + sum(is_internal_posterior(:)));

for left_side = [true false]
    for i=1:N_chordae

        left_side
        i

        errors = zeros(size(epsilon_vals)); 


        for ep_idx = 1:length(epsilon_vals)

            ep = epsilon_vals(ep_idx); 

            % make a new structure for the perturbation 
            valve_perturbation   = valve;  
            valve_perturbation.anterior.X  = valve.anterior.X  + ep*valve_Z.anterior.X; 
            valve_perturbation.posterior.X = valve.posterior.X + ep*valve_Z.posterior.X; 

            valve_perturbation.anterior.chordae.C_left  = valve.anterior.chordae.C_left  + ep*valve_Z.anterior.chordae.C_left; 
            valve_perturbation.anterior.chordae.C_right = valve.anterior.chordae.C_right + ep*valve_Z.anterior.chordae.C_right; 

            % eval the difference eqns on the perturbation 
            [F_anterior_perturbed F_posterior_perturbed F_chordae_left_perturbed F_chordae_right_perturbed] = difference_equations_bead_slip_attached(valve_perturbation); 

            F_perturbed_linearized = linearize_internal_points_bead_slip_attached(valve, F_anterior_perturbed, F_posterior_perturbed, F_chordae_left_perturbed, F_chordae_right_perturbed); 

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
 













