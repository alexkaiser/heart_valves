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


% 'solution immediately before exit'
% [X,S,T,N,du,dv,p_0] = unpack_params(params); 
% X(:,j,k)

if any(isnan(params.X(:,j,k)))
    error('STILL HAVE NANs on function exit, (j,k) = (%d,%d)\n', it, j, k);
end 

if any(isnan(params.X(:,j,k)))
    error('STILL HAVE INFs on function exit, (j,k) = (%d,%d)\n', it, j, k);
end 




