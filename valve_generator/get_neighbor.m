function [X_nbr R_nbr idx_chordae left] = get_neighbor(leaflet, j_nbr, k_nbr)
% 
% Returns the neighbor to a point on 
% By convention, if j==0 or k==0 then points off the leaflet are used 
% 
% leaflet.chordae must be nonempty 
% If nbr is 
% then the appropriate location on the chordae tree is returned 
% 

X = leaflet.X; 
R = leaflet.X; 

if leaflet.radial_and_circumferential
    error('Radial and circumferential fibers not implemented')
end 

idx_chordae = 0; 
left = false; 

if (j_nbr > 0) && (k_nbr > 0) 
    
    X_nbr = X(:,j_nbr,k_nbr); 
    R_nbr = R(:,j_nbr,k_nbr); 
    
elseif (j_nbr == 0) && (k_nbr == 0) 
    
    error('there is no zero zero direction neighbor'); 
     
else
    
    % chordae part 
    if (j_nbr == 0) && (k_nbr > 0)
        
        % use the left chordae 
        left = true; 
        [m max_internal] = size(leaflet.chordae.C_left); 
        
        % if the tree and the chordae overlapped 
        % what would the index (in tree numbering) be? 
        tree_idx_on_leaflet = max_internal + k_nbr; 
        
        % then take the parent index of that number in chordae variables 
        idx_chordae = floor(tree_idx_on_leaflet /2);  
        
        X_nbr = leaflet.chordae.C_left(:,idx_chordae); 
        R_nbr = leaflet.chordae.Ref_l(:,idx_chordae); 
        
    elseif (k_nbr == 0) && (j_nbr > 0)
        
        % use the right chordae 
        [m max_internal] = size(leaflet.chordae.C_right); 
        tree_idx_on_leaflet = max_internal + j_nbr; 
        idx_chordae = floor(tree_idx_on_leaflet/2); 
        X_nbr = leaflet.chordae.C_right(:,idx_chordae); 
        R_nbr = leaflet.chordae.Ref_r(:,idx_chordae); 
        
    else 
        error('should be impossible to get here'); 
    end 
end 





