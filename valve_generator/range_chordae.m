function range = range_chordae(chordae, idx_chordae, tree_idx)
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
 
if idx_chordae == 0
    range = [];
else 
    range = chordae(tree_idx).min_global_idx + 3*(idx_chordae-1) + (0:2);
end 

