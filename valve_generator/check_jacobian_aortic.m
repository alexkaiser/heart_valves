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

N = 16*3; 

% Initialize structures  
valve = initialize_valve_data_structures_aortic_generic(N); 

rand('twister',76599)

epsilon_vals = 10.^(-1:-1:-8); 

errors = zeros(size(epsilon_vals)); 

leaflet = valve.leaflets(1); 

% eval the difference eqns on the perturbation 
F  = leaflet.diff_eqns(leaflet); 

% jacobian does not change 
J = leaflet.jacobian(leaflet); 

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

Z_linearized = linearize_internal_points(leaflet_Z); 



fprintf('eps\t & Taylor series remainder\n'); 


for i = 1:length(epsilon_vals)
    
    ep = epsilon_vals(i); 
    
    % make a new structure for the perturbation 
    leaflet_perturbation   = leaflet;  
    leaflet_perturbation.X = leaflet.X + ep*leaflet_Z.X; 
    
    % eval the difference eqns on the perturbation 
    F_perturbed = leaflet.diff_eqns(leaflet_perturbation); 

    errors(i) = norm(F_perturbed - F - ep*J*Z_linearized, 2); 
    
    fprintf('%.1e\t & %e \\\\ \n \\hline \n', ep, errors(i)); 

end 

fprintf('\n\n\n\n'); 

fig = figure; 
loglog(epsilon_vals, errors, 'k-*'); 
hold on 
loglog(epsilon_vals, epsilon_vals.^2, 'k--'); 
axis([min(epsilon_vals) max(epsilon_vals) 1e-18, 1e8])
xlabel('\epsilon')
legend('|r(\epsilon)|', '\epsilon^2', 'Location', 'NorthWest')
printfig(fig, 'jacobian_convergence')


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
                
                % eval the difference eqns on the perturbation 
                F_perturbed = leaflet.diff_eqns(leaflet_perturbation); 

                diffs = F_perturbed - F - ep*J*Z_linearized; 

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

