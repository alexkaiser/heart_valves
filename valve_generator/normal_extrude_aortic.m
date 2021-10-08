function X_extruded = normal_extrude_aortic(leaflet, extrude_length)
    
    % takes the leaflet and extrudes normally 
    % 

    X           = leaflet.X; 
    j_max       = leaflet.j_max; 
    k_max       = leaflet.k_max; 
    is_internal = leaflet.is_internal; 
    is_bc       = leaflet.is_bc; 
    
    X_extruded  = zeros(size(X));  
    
%     if isfield(params, 'center_extrusion') && params.center_extrusion
%         
%         if params.num_copies ~= 3
%             error('cneter extrusion only implemted with 3 layers') 
%         end 
%         
%         extrude_length = (params.copy - 2) * params.ds_extrude; 
%     else 
%         extrude_length = (params.copy - 1) * params.ds_extrude; 
%     end 
    
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) || is_bc(j,k) 
                
                [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k); 
                    
                normal = cross(X(:,j_plus__1,k) - X(:,j_minus_1,k), X(:,j,k_plus__1) - X(:,j,k_minus_1));                     
                normal = normal / norm(normal); 
                
                if isfield(leaflet, 'extrusion_out') && leaflet.extrusion_out
                    normal = -normal; 
                end 
                
%                 if (j==j_max) && (k==k_max)
%                     'pause'
%                 end 
    
                if isfield(leaflet, 'fused_commissure') && leaflet.fused_commissure                     
                    if ~isfield(leaflet, 'fused_comm_idx')
                        error('must supply fused_comm_idx if leaflet.fused_commissure is true')
                    end 
                    
                    if leaflet.fused_comm_idx ~= 3
                        error('not implemented')
                    end 
                    
                    if (j==j_max) && (k==k_max) % comm point 
                        normal = [-1; 0; 0]; % just point straight in here 
                    end 
                    
                end 

                % if leaflet is being pinched
                % at comm point, ignore extrusion 
                % this point is removed anyway when thicken
                if isfield(leaflet, 'N_to_pinch') && (leaflet.N_to_pinch > 0)
                    if (mod(j, leaflet.N_each) == 0) && (k == k_max)
                        if any(isnan(normal))
                            normal = [0; 0; 0]; 
                        end 
                    end 
                end 
                
                X_extruded(:,j,k) = X(:,j,k) - normal * extrude_length; 
                
            end 
        end 
    end 
    
    if any(any(any(isnan(X_extruded))))
        error('Something not set in aortic valve extrustion'); 
    end 
   
end 
