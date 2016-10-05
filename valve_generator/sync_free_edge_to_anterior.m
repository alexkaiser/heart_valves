function X = sync_free_edge_to_anterior(anterior, posterior, X)
% 
% Synchronizes the free edge of the posterior leaflet to the anterior leaflet 
% 
% 

free_edge_idx_left      = anterior.free_edge_idx_left; 
free_edge_idx_right     = anterior.free_edge_idx_right; 

is_bc = posterior.is_bc; 

for i=1:size(free_edge_idx_left, 1)
    j = free_edge_idx_left(i,1); 
    k = free_edge_idx_left(i,2); 
    X(:,j,k) = anterior.X(:,j,k);  
    
    if ~is_bc(j,k)
        error('Synchronizing to location that is not a boundary condition point')
    end 
    
end

for i=1:size(free_edge_idx_right, 1)
    j = free_edge_idx_right(i,1); 
    k = free_edge_idx_right(i,2); 
    X(:,j,k) = anterior.X(:,j,k);  
    
    if ~is_bc(j,k)
        error('Synchronizing to location that is not a boundary condition point')
    end 
    
end 

