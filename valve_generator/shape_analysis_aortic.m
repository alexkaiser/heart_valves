function [Eh, Gh, min_height_center, orifice_area, free_edge_length] = shape_analysis_aortic(leaflet)

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 
X      = leaflet.X; 

R_u    = leaflet.R_u; 
R_v    = leaflet.R_v; 

if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    N_leaflets = 3; 
end 

Eh = zeros(N_leaflets,1); 
Gh = zeros(N_leaflets,1); 
min_height_center = zeros(N_leaflets,1); 
orifice_area = 0.0; 
free_edge_length = zeros(N_leaflets,1);


% height 
for comm_idx = 1:N_leaflets

    min_idx = (comm_idx-1)*N_each;         
   
    center_idx_j = min_idx + N_each/2

    j = center_idx_j; 
    
    for k=2:k_max
    
        j_nbr_tmp = j;
        k_nbr_tmp = k-1; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = X(:,j,k);
        X_nbr = X(:,j_nbr,k_nbr); 

        Gh(comm_idx) = Gh(comm_idx) + norm(X_temp - X_nbr);        
    end 
    
    Eh(comm_idx) = X(3,j,k_max); 
    
    min_height_center(comm_idx) = min(X(3,j,:));
    
end



% free edge lengths 
for comm_idx = 1:N_leaflets

    min_idx = (comm_idx-1)*N_each;         
    
    for j=(1 + min_idx):(N_each + min_idx)
        k = k_max; 

        j_nbr_tmp = j-1; 
        k_nbr_tmp = k; 
        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end

        X_temp = X(:,j,k);
        X_nbr = X(:,j_nbr,k_nbr); 

        free_edge_length(comm_idx) = free_edge_length(comm_idx) + norm(X_temp - X_nbr);        
        
    end
    
end


orifice_area = polyarea(X(1,:,k_max), X(2,:,k_max)); 












