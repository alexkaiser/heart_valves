function range = range_chordae(total_internal, N_chordae, idx_chordae, left)
% 
% Returns the range of indices corresponding to 
% 
% Input
%     total_internal    Total number of scalar points in the leaflet
%                       Does not include b.c.
%     N_chordae         Number of (vector) points in the tree of chordae 
%     idx_chordae       Current index in the (vector index) direction in the chordae 
%     left              Whether to use the left side of the chordae 
% 
% Output
%     range             Length three array of global indices 
% 
 
if left
    range = total_internal + 3*(idx_chordae-1) + (1:3);  
else
    range = total_internal + 3*N_chordae + 3*(idx_chordae-1) + (1:3);  
end 

% if the index requested is zero, then it is a papillary muscle 
% which is a b.c. 
if idx_chordae == 0
    range = []; 
end 












