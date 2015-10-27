function X = build_reference_surface(filter_params)
    %
    % Builds reference coffee cone surface 
    % 
    % Mesh has a triangular layout 
    % j==1 and k==1 are assumed to be the free edge 
    % j+k == N+2 is the valve ring 
    % 
    % 
    % 

    N = filter_params.N; 
    r = filter_params.r; 
    h = filter_params.h; 

    X      = zeros(3,N+1,N+1); 
    X_flat = zeros(2,N+1,N+1); 

    mesh = linspace(-pi/2,pi/2,N+1); 
    ring_half = [r*cos(mesh); r*sin(mesh); h*ones(size(mesh))]; 

    % set the valve ring 
    for j=1:N+1
        k = (N+2) - j; % ring coodinates 
        X_flat(:,j,k) = cone_filter_inv(ring_half(:,j), filter_params); 
        X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
    end 

    % fill in the 3d array 
    for j=1:N
        for k=1:N
            
            % in the triangle? 
            if ((j+k) < (N+2))
                X_flat(:,j,k) = compute_intersection(X_flat, j, k, filter_params); 
                X(:,j,k)      = cone_filter(X_flat(1,j,k), X_flat(2,j,k), filter_params); 
            end

        end 
    end

end 




function coords = compute_intersection(X_flat, j, k, filter_params)
% 
% Computes the intersection of the rays which have coordinates j,k
% 
% Input: 
%     X_flat          2d surface including valve ring 
%     j,k             Indices to fill 
%     filter_params   paramters for the filter 
% 
% Output: 
%     coords          2d vector for intersection 
% 

    a = filter_params.a; 
    
    
    % find the indices corresponding diagonal points which lie on the ray
    j_left  = filter_params.N + 2 - k; 
    k_right = filter_params.N + 2 - j; 
    
    % intersection is a solution to a linear system 
    A = [a + X_flat(1,j_left,k), a - X_flat(1,j,k_right); 
             X_flat(2,j_left,k),   - X_flat(2,j,k_right)]; 

    rhs = [X_flat(1,j_left,k) - X_flat(1,j,k_right); 
           X_flat(2,j_left,k) - X_flat(2,j,k_right); ]; 
    
    soln = A \ rhs; 
    s = soln(1); 
    
    coords = [-a;0]*s + (1 - s)*X_flat(:,j_left,k); 

    t = soln(2); 
    coords_alt = [a;0]*t + (1 - t)*X_flat(:,j,k_right);
    
    tol = 1e2 * eps; 
    if norm(coords - coords_alt) > tol 
        error('two formulas from line disagree'); 
    end 
    
end 




















