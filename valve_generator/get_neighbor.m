function [X_nbr R_nbr] = get_neighbor(leaflet, j, k, j_nbr, k_nbr)
% 
% Returns the neighbor to a point on the leaflet 
% 
% If j_nbr,k_nbr is on the leaflet or a b.c. then its coordinates are returned 
% 
% If not, then it is checked whether j,k has chordae tendineae, 
% and chordae are returned 
% 

X = leaflet.X; 
R = leaflet.R; 

if leaflet.radial_and_circumferential
    error('Radial and circumferential fibers not implemented')
end 

% idx_chordae = 0; 
% left = false; 

if (j_nbr > 0) && (k_nbr > 0) && (leaflet.is_internal(j_nbr,k_nbr) || leaflet.is_bc(j_nbr,k_nbr))
    
    X_nbr = X(:,j_nbr,k_nbr); 
    R_nbr = R(:,j_nbr,k_nbr); 
    
% elseif leaflet.chordae_idx_left(j,k)
%     
%     [m max_internal] = size(leaflet.chordae.C_left); 
%     
%     % if the leaves of the tree were stored, tree index would have this value 
%     leaf_idx = leaflet.chordae_idx_left(j,k) + max_internal; 
%     
%     % then take the parent index of that number in chordae indexing  
%     idx_chordae = floor(leaf_idx /2);  
% 
%     X_nbr = leaflet.chordae.C_left(:,idx_chordae); 
%     R_nbr = leaflet.chordae.Ref_l (:,idx_chordae);     
%     
% elseif leaflet.chordae_idx_right(j,k)
%     
%     [m max_internal] = size(leaflet.chordae.C_right); 
%     leaf_idx = leaflet.chordae_idx_right(j,k) + max_internal; 
%     idx_chordae = floor(leaf_idx /2);  
%     X_nbr = leaflet.chordae.C_right(:,idx_chordae); 
%     R_nbr = leaflet.chordae.Ref_r  (:,idx_chordae);  
    
else 
    error('No valid neighbor found')
end 
    
    
    
% elseif (j_nbr == 0) && (k_nbr == 0) 
%     
%     error('there is no zero zero direction neighbor'); 
%      
% else
%     
%     % chordae part 
%     if (j_nbr == 0) && (k_nbr > 0)
%         
%         % use the left chordae 
%         left = true; 
%         [m max_internal] = size(leaflet.chordae.C_left); 
%         
%         % if the tree and the chordae overlapped 
%         % what would the index (in tree numbering) be? 
%         tree_idx_on_leaflet = max_internal + k_nbr; 
%         
%         % then take the parent index of that number in chordae variables 
%         idx_chordae = floor(tree_idx_on_leaflet /2);  
%         
%         X_nbr = leaflet.chordae.C_left(:,idx_chordae); 
%         R_nbr = leaflet.chordae.Ref_l(:,idx_chordae); 
%         
%     elseif (k_nbr == 0) && (j_nbr > 0)
%         
%         % use the right chordae 
%         [m max_internal] = size(leaflet.chordae.C_right); 
%         tree_idx_on_leaflet = max_internal + j_nbr; 
%         idx_chordae = floor(tree_idx_on_leaflet/2); 
%         X_nbr = leaflet.chordae.C_right(:,idx_chordae); 
%         R_nbr = leaflet.chordae.Ref_r(:,idx_chordae); 
%         
%     else 
%         error('should be impossible to get here'); 
%     end 
% end 





