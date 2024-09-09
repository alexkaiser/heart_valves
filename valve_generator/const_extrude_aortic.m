function X_extruded = const_extrude_aortic(leaflet, extrude_length)
    
    % takes the leaflet and extrudes normally 
    % 

    X           = leaflet.X; 
    j_max       = leaflet.j_max; 
    k_max       = leaflet.k_max; 
    is_internal = leaflet.is_internal; 
    is_bc       = leaflet.is_bc; 
    N           = leaflet.N; 
    N_each      = leaflet.N_each; 
    N_leaflets  = leaflet.N_leaflets; 
    
    X_extruded  = nan * zeros(size(X));  
    
    tol = 1e-10; 
            
    idx_comm_below = 1;  
    idx_comm_above = j_max; 


    chord = X(:,idx_comm_above,1) - X(:,idx_comm_below,1); 

    if abs(chord(3)) > tol 
        error('chord not horizontal'); 
    end 

    % zero it for good measure 
    chord(3) = 0; 
    chord = chord / norm(chord);

    % orthogonal to given chord         
    normal = [-chord(2); chord(1); 0]; 

    if isfield(leaflet, 'extrusion_out') && leaflet.extrusion_out
        normal = -normal; 
    end 

    for j = 1:j_max
        for k=1:k_max
            if is_internal(j,k) || is_bc(j,k) 
                X_extruded(:,j,k) = X(:,j,k) + normal * extrude_length; 
            end 
        end 
    end 
                 
    if any(any(any(isnan(X_extruded))))
        error('Something not set in aortic valve extrustion'); 
    end 
   
end 
