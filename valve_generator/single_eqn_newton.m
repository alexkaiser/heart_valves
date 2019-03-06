function [params pass] = single_eqn_newton(params,filter_params,j,k,tol,max_it)
%
% Newton solve on a single j,k set of equations 
% 
% If fails to converge, X is not changed and the intial values are returned 
% This is potentially questionable... 
% 
% params,
% j
% k
% tol
% max_it
% jacobian

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


X_initial = params.X(:,j,k); 

F = difference_equations_component(params, filter_params,j,k); 
residual = norm(F); 

pass = true; 
it = 0; 


while residual > tol 
    
    it = it + 1; 
    
    J = build_jacobian_component(params, filter_params,j,k); 
    
    % solve the system, just use backslash it's 3 by 3
    soln = J \ (-F); 
    
    % add in to get the next iterate 
    params.X(:,j,k) = params.X(:,j,k) + soln(1:3); 
    
    F = difference_equations_component(params, filter_params,j,k); 
    residual = norm(F); 
    
    % fprintf('it = %d, \tresidual = %e\n', it, residual); 
    
    if it > max_it 
        warning('Newton solve failed to converge in %d iterations, (j,k) = (%d,%d)\n', it, j, k);
        pass = false; 
        break; 
    end 
    
    if any(isnan(params.X(:,j,k)))
        warning('NAN problems: Newton solve failed to converge in %d iterations, (j,k) = (%d,%d)\n', it, j, k);
        pass = false; 
        break; 
    end 
    
    if any(isinf(params.X(:,j,k))) 
        warning('INF problems: Newton solve failed to converge in %d iterations, (j,k) = (%d,%d)\n', it, j, k);
        pass = false; 
        break; 
    end
    
end 



% current version will 

if pass 
    if (j<3) && (k<3)
        fprintf(1, 'Newton solve converged in %d iterations with residual = %g, (j,k) = (%d,%d)\n', it, residual,j,k);
    end 
% elseif (residual < (0.1*initial_residual)) && (~any(isnan(X(:,j,k))) && (~any(isinf(X(:,j,k))))
%     warning('Despite failure, residual reduced and current iterate used.');
else 
    params.X(:,j,k) = X_initial; 
end

if any(isnan(params.X(:,j,k)))
    error('STILL HAVE NANs on function exit, (j,k) = (%d,%d)\n', it, j, k);
end 

if any(isnan(params.X(:,j,k)))
    error('STILL HAVE INFs on function exit, (j,k) = (%d,%d)\n', it, j, k);
end 




