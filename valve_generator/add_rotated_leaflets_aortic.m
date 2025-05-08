function valve = add_rotated_leaflets_aortic(valve)

N_leaflets = valve.leaflets(1).N_leaflets; 

if length(valve.leaflets) ~= 1 
    error('rotation implmented only with one initial leaflet'); 
end 

if N_leaflets == 1
    warning('No rotation of one single leaflet')
    return; 
end 

th = 2*pi/N_leaflets; 

k_max = valve.leaflets(1).k_max; 

for n = 2:N_leaflets
    
    valve.leaflets(n) = valve.leaflets(1);
    
    % rotate all cols 
    for k = 1:k_max
        valve.leaflets(n).X(:,:,k) = rotation_matrix_z((n-1)*th) * valve.leaflets(n).X(:,:,k);
    end 
    
end 









