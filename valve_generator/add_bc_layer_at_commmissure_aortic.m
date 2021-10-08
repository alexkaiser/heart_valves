function leaflet = add_bc_layer_at_commmissure_aortic(leaflet)
    % labels the radial fiber to be boundary conditions 
    % adjacent in the j direction of existing boundary conditions  
    % this means that the adjacent radial fiber will be assigned to be a target 
    % a single circumferential link on each fiber will be removed 
    % and that link will become part of the commissure

    j_max            = leaflet.j_max; 
    k_max            = leaflet.k_max; 
    is_bc            = leaflet.is_bc; 
    is_internal_temp = leaflet.is_internal; 
    is_bc_temp       = is_bc; 
    
    for j=1:j_max
        for k=2:k_max
            if is_bc(j,k) 
                [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k); 

                is_bc_temp(j_minus_1,k)       = true; 
                is_bc_temp(j_plus__1,k)       = true; 
                is_internal_temp(j_minus_1,k) = false; 
                is_internal_temp(j_plus__1,k) = false; 
            end 
        end 
    end 

    leaflet.is_bc       = is_bc_temp; 
    leaflet.is_internal = is_internal_temp; 

end 
