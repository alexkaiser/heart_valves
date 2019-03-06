function [] = check_single_fiber()

    % 
    % hacked together test of tensions and jacobians 
    % with slip models 
    % 
    % makes a single fiber, computes the 

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

    N = 8; 

    X = zeros(3,N+1); 

    % boundary condition points 
    X_current(:,1) = zeros(3,1); 
    X_current(:,N+1) = ones(3,1); 

    % points 2 is a "free edge" points and gets a linear constitutive model 
    % other points in 


    k_max = N; 

    for i=2:N
        X_current(:,i,1) = (1 - i/(N+1)) * X_current(:,1,1) + (i/(N+1)) * X_current(:,N+1,1); 
    end 

    ref_frac = 0.7; 
    Ref = NaN * zeros(size(X)); 

    % first three points have a resting location 
    Ref(:,1:3) = X_current(:,1:3); 

    kappa = 1; 
    du = 1/(N+1); 

    % t type first 
    [F, T] = diff_eqns_one_sided(X_current, Ref, kappa, ref_frac, du, N); 

    J = zeros(3*(N+1),3*(N+1)); 

    % linear tension Jacobians 
    j=2; 
    for j_nbr = [j-1,j+1]
        X = X_current(:,j); 
        R = Ref(:,j); 
        X_nbr = X_current(:,j_nbr); 
        R_nbr = Ref(:,j_nbr); 

        J_tension = tension_linear_tangent_jacobian(X, X_nbr, R, R_nbr, kappa, ref_frac);

        range_current = 3*(j-1) + (1:3);     
        range_nbr = 3*(j_nbr-1) + (1:3); 

        J(range_current, range_current) = J(range_current, range_current) + J_tension;

        if (j_nbr > 1) && (j_nbr < (N+1))
            J(range_current, range_nbr) = J(range_current, range_nbr) + -J_tension; 
        end 

    end

%     % slip model jacobians 
%     for j=3:N
%         for j_nbr = [j-1,j+1]
%             X = X_current(:,j); 
%             X_nbr = X_current(:,j_nbr);
% 
%             tension = (1/du) * 0;  % * T.val; 
%             grad_tension  = 0*(1/du) * T.G; 
%             tangent = (X_nbr - X) / norm(X_nbr - X); 
%             J_tangent = tangent_jacobian(X, X_nbr);
% 
%             range_current = 3*(j-1) + (1:3); 
% 
%             range_nbr = 3*(j_nbr-1) + (1:3); 
% 
%             J(range_current, range_current) = J(range_current, range_current) + tension * J_tangent;
% 
%             if (j_nbr > 1) && (j_nbr < (N+1))
%                 J(range_current, range_nbr) = J(range_current, range_nbr) + -tension * J_tangent; 
%             end 
% 
%             j_edge = T.j; 
%             Jac = grad_tension * tangent';
%             range_edge = 3*(j_edge-1) + (1:3); 
%             if (j_edge > 1) && (j_edge < (N+1))
%                 J(range_current, range_edge) = J(range_current, range_edge) + Jac; 
%             end 
% 
%             j_tension = T.j_nbr; 
%             range_tension = 3*(j_tension-1) + (1:3); 
%             if (j_tension > 1) && (j_tension < (N+1))
%                 J(range_current, range_tension) = J(range_current, range_tension) + -Jac;         
%             end 
%         end 
%     end 



    Z = rand(size(X_current)); 
    Z(:,1)   = zeros(3,1); 
    Z(:,N+1) = zeros(3,1); 

    epsilon_vals = 10.^(-1:-1:-8); 
    errors = zeros(size(epsilon_vals)); 

    for i = 1:length(epsilon_vals)


        ep = epsilon_vals(i); 

        P = X_current + ep*Z; 

        F_perturbed = diff_eqns_one_sided(P, Ref, kappa, ref_frac, du, N); 

        errors(i) = norm(F_perturbed(:) - F(:) - ep*J*Z(:), 2); 

        fprintf('%e\t | %e \n', ep, errors(i)); 

    end 

end 

function [F T] = diff_eqns_one_sided(X_current, Ref, kappa, ref_frac, du, N)

    F = zeros(size(X_current)); 

    % first index is connected to b.c. with a linear spring 
    j=2; 
    for j_nbr = [j-1,j+1]
        X = X_current(:,j); 
        R = Ref(:,j); 
        X_nbr = X_current(:,j_nbr); 
        R_nbr = Ref(:,j_nbr); 
        F(:,j) = F(:,j) + tension_linear(X, X_nbr, R, R_nbr, kappa, ref_frac); 
    end 

    j=2; 
    j_nbr = 3; 

    X = X_current(:,j); 
    R = Ref(:,j); 
    X_nbr = X_current(:,j_nbr); 
    R_nbr = Ref(:,j_nbr);

    k = []; 
    k_nbr = []; 

    T = get_linear_tension_struct(X, X_nbr, R, R_nbr, kappa/du, ref_frac, j, k, j_nbr, k_nbr); 

    for j=3:N    
        for j_nbr = [j-1,j+1]
            X = X_current(:,j); 
            X_nbr = X_current(:,j_nbr);
            % F(:,j) = F(:,j) + (1/du)*T.val * (X_nbr-X)/norm(X_nbr-X); 
            F(:,j) = F(:,j) + (1/du) *0* (X_nbr-X)/norm(X_nbr-X); 
        end 
    end 


end 





































