function X = cone_filter(xi, eta, filter_params)
%
% Evaluates the cone filter at (xi, eta)
% 
% Input 
%     xi, eta         Euclidean coordinates in plane 
%     filter_params   Struct with cone filter parameters 
%     
% Output
%     X               Vector of components 
% 


    a = filter_params.a; 
    r = filter_params.r; 
    h = filter_params.h; 

    tol = 1e2 * eps;

    % If on the pole return it
    if abs(eta) < tol && abs(xi-a) < tol 
        X = [0; a;0]; 
        return; 
    end 

    % If on the other pole return it
    if abs(eta) < tol && abs(xi+a) < tol 
        X = [0;-a;0]; 
        return; 
    end 

    % maximum angle in both coordinate systems 
    theta_0 = acos(-a / sqrt(a^2 + r^2 + h^2));

    % minimum angle in both coordinate systems to be on the front sheet
    integrand = @(phi) - sqrt( (1 - (a/r) * sin(phi)).^2 + (h/r)^2) ./ (1 - (2*a/r) * sin(phi) + (a/r)^2 + (h/r)^2) ; 
    theta_min = theta_0 + quadgk(integrand,0,pi/2,'RelTol',tol,'AbsTol',tol);
    
    % right polar coordinate angle first 
    theta_right = atan2(eta, xi - a);

    % Other angle 
    % use reflected coordinates, (-xi,eta)
    % and the (a,0) polar coordinate origin 
    theta_left = atan2(eta, -xi - a);

    % are we in the triangle portion, which simply gets rotated out? 
    if (theta_0 <= theta_left) && (theta_0 <= theta_right)

        % we are on the triangle, rotate out 
        psi = asin(r / sqrt(h^2 + r^2)); 

        X = zeros(3,1); 
        X(1) = sin(psi) * eta; 
        X(2) =             xi; 
        X(3) = cos(psi) * eta;     


    % If angle is less than the leftmost ray on the front 
    % then point is on the right, front sheet 
    elseif (theta_min <= theta_right) && (theta_right <= theta_0)
        X = cone_filter_right_cone(xi, eta, filter_params); 
        return 

    % are we on the back of the right sheet 
    elseif ((2*theta_min - theta_0) <= theta_right) && (theta_right <= theta_min)

        % translate, rotate and reflect to use the other coordinates 
        
        % reset origin 
        xi = xi - a; 

        % rotate axis to positive eta axis
        val = rotation_matrix(theta_min - pi/2) * [xi; eta]; 
        xi  = val(1); 
        eta = val(2); 
        
        % reflect over eta axis, negate xi 
        xi = -xi; 

        % rotate back 
        val = rotation_matrix(-(theta_min - pi/2)) * [xi; eta]; 
        xi  = val(1); 
        eta = val(2); 
        
        % unset origin
        xi = xi + a; 
        
        % apply filter map 
        X = cone_filter_right_cone(xi, eta, filter_params); 
        
        % negate x for the back sheet 
        X(1) = -X(1); 
        
        return 


    % left front sheet 
    elseif (theta_min <= theta_left) && (theta_left <= theta_0)
        % apply the reflection in xi to get the corresponding preimage 
        X = cone_filter_right_cone(-xi, eta, filter_params); 

        % reflect back in y in the 3d setting  
        X(2) = -X(2); 
        return; 

    elseif ((2*theta_min - theta_0) <= theta_left) && (theta_left <= theta_min)
        
        % same as above but various signs are swapped 
        
        % reset origin 
        xi = xi + a; 

        % rotate axis to positive eta axis
        val = rotation_matrix(-(theta_min - pi/2)) * [xi; eta]; 
        xi  = val(1); 
        eta = val(2); 
        
        % reflect over eta axis, negate xi 
        xi = -xi; 

        % rotate back 
        val = rotation_matrix(theta_min - pi/2) * [xi; eta]; 
        xi  = val(1); 
        eta = val(2); 
        
        % unset origin
        xi = xi - a; 
        
        % apply filter map using the LEFT map 
        X = cone_filter_right_cone(-xi, eta, filter_params); 
        
        % negate in y for the left side 
        X(2) = -X(2); 
        
        % negate x for the back sheet 
        X(1) = -X(1); 
        
        return 
        
    else 
        theta_min
        theta_0
        'min back sheet boundary'
        2*theta_min - theta_0
        theta_right
        theta_left
        error('should not have gotten here, coordinates not in any acceptable regions'); 
    end 

end 







