function leaflet = internal_points_to_2d(v_linearized, leaflet)
%
%  Takes the internal values in X which are arranged in linear order
%  And places them back in the 3d vector array in leaflet
%  

idx = 1; 
N = leaflet.N; 

if leaflet.radial_and_circumferential
    error('Radial and circumferential fibers not implemented')
end 

% here k is required to be the outer loop 
for k=1:N
    for j=1:N
        if leaflet.is_internal(j,k)
            leaflet.X(:,j,k) = v_linearized(idx + (0:2)); 
            idx = idx + 3; 
        end 
    end 
end 


% copy chordae if length allows 
if length(v_linearized)> idx
    
    C_left   = leaflet.chordae.C_left; 
    C_right  = leaflet.chordae.C_right; 

    [m N_chordae] = size(C_left); 
    total_internal = 3*N*(N+1)/2; 

    idx = total_internal + 1; 
    for i=1:N_chordae
        C_left(:,i)  = v_linearized(idx + (0:2));  
        idx = idx + 3; 
    end 

    idx = total_internal + 3*N_chordae + 1; 
    for i=1:N_chordae
        C_right(:,i) = v_linearized(idx + (0:2));  
        idx = idx + 3; 
    end 

    leaflet.chordae.C_left  = C_left; 
    leaflet.chordae.C_right = C_right; 

end 
