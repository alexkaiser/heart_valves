function order_check(path, N_values)
% 
% 
% 
% 



iterations = length(N_values) - 1; 

diff_1    = zeros(iterations,1); 
diff_2    = zeros(iterations,1); 
diff_inf  = zeros(iterations,1); 


for i = 1:iterations 

    N_coarse = N_values(i); 
    data_coarse = sprintf('%s/mitral_tree_%d_final_data', path, N_coarse); 
    load(data_coarse) 

    valve_coarse       = valve; 
    leaflet_coarse     = valve_coarse.leaflets(1); 
    X_coarse           = leaflet_coarse.X; 
    j_max_coarse       = leaflet_coarse.j_max; 
    k_max_coarse       = leaflet_coarse.k_max; 
    du_coarse          = leaflet_coarse.du; 
    is_internal_coarse = leaflet_coarse.is_internal;
    
    N_fine = 2*N_coarse; 
    data_fine = sprintf('%s/mitral_tree_%d_final_data', path, N_fine); 
    load (data_fine) 
    
    
    valve_fine       = valve; 
    leaflet_fine     = valve_fine.leaflets(1); 
    X_fine           = leaflet_fine.X; 
    du_fine          = leaflet_fine.du; 
    j_max_fine       = leaflet_fine.j_max; 
    k_max_fine       = leaflet_fine.k_max; 
    is_internal_fine = leaflet_fine.is_internal; 
    is_bc_fine       = leaflet_fine.is_bc;
    
    
    
    
    points_written_coarse = 0; 
    points_written_fine   = 0; 
    
    % points on X_coarse that have four valid neighbors 
    % on fine mesh for comparison 
    X_coarse_valid = zeros(size(X_coarse)); 
    
    % fine mesh interpolated to align with coarse mesh 
    X_fine_interp  = zeros(size(X_coarse)); 
    
    % for all internal points on coarse mesh 
    % find coordinates in u,v plane 
    % interpolate fine mesh to this location 
    % with bilinear interpolation 
    % only do so if all four points on the fine mesh are valid 
    for j=1:j_max_coarse
        for k=1:k_max_coarse
            if is_internal_coarse(j,k)
    
                % values of u,v at which function is defined 
                u = du_coarse * j; 
                v = du_coarse * k; 
    
                j_fine_below = floor(u/du_fine); 
                k_fine_below = floor(v/du_fine); 
                
                % make sure that we are looking at a valid index 
                if (0 < j_fine_below) && (j_fine_below <= j_max_fine) && ... 
                   (0 < k_fine_below) && (k_fine_below <= k_max_fine) && ... 
                   (is_internal_fine(j_fine_below, k_fine_below) || is_bc_fine(j_fine_below, k_fine_below)) 
               
                    j_nbr_tmp = j_fine_below + 1; 
                    k_nbr_tmp = k_fine_below + 1; 
               
                    % check if have valid neighbor
                    % and get indices reduced for periodicity etc  
                    [valid j_fine_above k_fine_above] = get_indices(leaflet_fine, j_fine_below, k_fine_below, j_nbr_tmp, k_nbr_tmp);  
                    
                    if valid 
                        
                        % have all necessary data to compare norms at this point 
                        % copy data into coarse 
                        X_coarse_valid(:,j,k) = X_coarse(:,j,k); 
                        points_written_coarse = points_written_coarse + 1; 
                        
                        if (any(isnan(X_coarse_valid(:,j,k))) || any(X_coarse_valid(:,j,k) == [0;0;0]))
                            error('nan or zero problems in coarse')
                        end 
                        
                        % compute fractions of distances, which are the interpolation coefficients 
                        % this is whatever remains after taking floor 
                        
                        % temp weight to only use bottom 
                        % this should still congeverge just with no particular order 
                        
                        
                        s = 1; %(u/du_fine) - j_fine_below; 
                        t = 1; %(v/du_fine) - k_fine_below; 
                        
                        X_fine_interp(:,j,k) = X_fine(:,j_fine_below,k_fine_below); 
                        
%                         X_fine_interp(:,j,k) = t     * (s * X_fine(:,j_fine_below,k_fine_below) + (1-s)*X_fine(:,j_fine_above,k_fine_below)) + ... 
%                                                (1-t) * (s * X_fine(:,j_fine_below,k_fine_above) + (1-s)*X_fine(:,j_fine_above,k_fine_above)); 
                        
                                           
                        if (any(isnan(X_fine_interp(:,j,k))) || any(X_fine_interp(:,j,k) == [0;0;0]))
                            error('nan or zero problems in fine')
                        end 
                                           
                        points_written_fine = points_written_fine + 1;                    
                    end 
                end 
                
            end 
        end 
    end 
    
    % very simple sanity check 
    if ~all(find(X_coarse_valid(:)) == find(X_fine_interp(:)))
        warning('nonzero indices do not match, check results')
    end 
    
    
    % For results to be comparable across vectors of different sizes
    % Need to scale norms so that they are double integrals 
    
    differences = X_coarse_valid(:) - X_fine_interp(:);  
    
    % du^2 measure for double integral 
    diff_1(i)    = du_coarse*du_coarse * norm(differences, 1); 
    diff_2(i)    = du_coarse           * norm(differences, 2); 
    diff_inf(i)  =                       norm(differences, 'inf'); 

end 

diff_1 
diff_2 
diff_inf 


fprintf('Comparisons of N,2N and 2N,4N points\n'); 
N = N_values(1); 
for i = 1:iterations-1

    order_1    = diff_1(i) / diff_1(i+1); 
    order_2    = diff_2(i) / diff_2(i+1);  
    order_inf  = diff_inf(i) / diff_inf(i+1);  

    fprintf('%d & %f & %f  &  %f    \\\\  \n \\hline \n ',  N,  order_1, order_2, order_inf);


    N = 2*N;
end















