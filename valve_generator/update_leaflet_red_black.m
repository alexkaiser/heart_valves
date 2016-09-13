function [params] = update_leaflet_red_black(params, filter_params, tol_local, max_it_local)
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
 

% Loop on internal in zero based indexing 
% Start at one, end at N

% Note that internal points are 1...N
% Zero is free edge 
% j+k = N+1 is valve ring 

for j=1:params.N
    for k=1:params.N
        if mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (params.N+2)
                params = single_eqn_newton(params,filter_params,j,k,tol_local,max_it_local); 
            end 
        end 
    end 
end

for j=1:params.N
    for k=1:params.N
        if ~mod(j+k,2)
            
            % in the triangle?
            if (j+k) < (params.N+2)
                params = single_eqn_newton(params,filter_params,j,k,tol_local,max_it_local); 
            end
            
        end 
    end 
end


    