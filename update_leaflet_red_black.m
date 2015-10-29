function [params] = update_leaflet_red_black(params, fiber_params, tol_local, max_it_local)
%
% Full leaflet and boundary sweep. 
% Solves the nonlinear difference equations at each component. 
% Updates params in red black ordering 
% 
% Input 
%     params             Current solution parameters 
%     tol_local          Error tolerance for local newton solve 
%     max_it_local       Max iterations on each component of the Newton solve 
%     triangle_domain    Only evalutates internal points if they are below the 
%                            line y + x = 1
% 
% Output 
%     params             Current solution parameters 
% 

% [X,alpha,beta,N,p_0,R] = unpack_params(params)
% 

% Loop on internal in zero based indexing 
% Start at one, end at N

% Note that internal points are 1...N
% Zero is free edge 
% j+k = N+1 is valve ring 

for j=1:N
    for k=1:N
        if mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (N+2)
                % params = single_eqn_newton(params,j+1,k+1,tol_local,max_it_local,jacobian,diff_eqns);
            end 
        end 
    end 
end

for j=0:N
    for k=0:N
        if ~mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (N+2)
                % params = single_eqn_newton(params,j+1,k+1,tol_local,max_it_local,jacobian,diff_eqns);
            end
            
        end 
    end 
end


    