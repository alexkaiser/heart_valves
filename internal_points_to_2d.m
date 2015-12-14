function params = internal_points_to_2d(v_linearized, params)
%
%  Takes the internal values in X which are arranged in linear order
%  And places them back in the 3d vector array in params
%  

idx = 1; 
N = params.N; 


% here k is required to be the outer loop 
for k=1:N
    for j=1:N
        % in the triangle?
        if (j+k) < (N+2)
            params.X(:,j,k) = v_linearized(idx + (0:2)); 
            idx = idx + 3; 
        end 
    end 
end 


% copy chordae 
if isfield(params, 'chordae') && ~isempty(chordae)
    
    [m N_chordae] = size(chordae.C_left); 
    total_internal = 3*N*(N+1)/2; 
    
    idx = total_internal + 1; 
    for i=1:N_chordae
        params.chordae.C_left(:,i)  = v_linearized(idx + (0:2));  
        idx = idx + 3; 
    end 
    
    idx = total_internal + N_chordae + 1; 
    for i=1:N_chordae
        params.chordae.C_right(:,i) = v_linearized(idx + (0:2));  
        idx = idx + 3; 
    end 
        
end 


