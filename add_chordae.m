function params = add_chordae(params, filter_params, k_0, k_multiplier, tree_frac, arbitrary_papillary_points)
% 
% Adds a chordae data structure to the current parameters 
% 
% Input: 
%     params           Current leaflet data structures 
%     filter_params    Cone filter parameters 
%     k_0              Spring constant from leaflet to first internal tree layer 
%     k_multiplier     Each level of the tree towards the
% 
% Output: 
%     params           Current leaflet data structures with chordae added 
% 


if arbitrary_papillary_points
    
    params.chordae = build_3d_trees_arbitrary_coords(params, filter_params, tree_frac, k_0, k_multiplier); 
    
else 
    left_papillary_flat  = [-filter_params.a; 0]; 
    right_papillary_flat = [ filter_params.a; 0]; 

    chordae_flat = build_2d_trees(params, left_papillary_flat, right_papillary_flat, tree_frac, filter_params); 
    chordae = build_3d_trees(chordae_flat, filter_params, k_0, k_multiplier); 

    params.chordae = chordae; 
end 
