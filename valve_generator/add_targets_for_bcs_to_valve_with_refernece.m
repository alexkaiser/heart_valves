function valve = add_targets_for_bcs_to_valve_with_refernece(valve, adjustment_length)
% 
% If not already included, adds necessary data structures to treat current valve configuration 
% as including target points for boundary conditions 
% 

debug = true; 

if isfield(valve, 'targets_for_bcs') && valve.targets_for_bcs
    return; 
end 
valve.targets_for_bcs = true; 

if ~exist('adjustment_length', 'var')
    adjustment_length = 1e-5; 
end 

leaflet = valve.leaflets(1); 

chordae_copy = leaflet.chordae; 

leaflet.targets_for_bcs = true; 

k_max = leaflet.k_max; 

% copy the bc values 
leaflet.X(:,:,k_max+1) = leaflet.X(:,:,k_max); 

% valve ring adjusted down in z component 
leaflet.X(3,:,k_max) = leaflet.X(3,:,k_max) - adjustment_length; 

% we've added a point here 
leaflet.k_max = k_max + 1; 
leaflet.ring_k_idx = leaflet.ring_k_idx + 1; 

% redo util arrays 
leaflet = get_util_arrays_bead_slip(leaflet, valve); 

% and update associated constants 

% Total number of internal leaflet coordinates (three times number of vertices)
leaflet.total_internal_leaflet    = 3*sum(leaflet.is_internal(:)); 

% Running total number of coordinates including trees 
% Updated as trees are added 
leaflet.total_internal_with_trees = 3*sum(leaflet.is_internal(:)); 

% spring constants and rest lengths do not need updating here 
% because they should never be read anyway 

% add chordae terms 
for tree_idx = 1:leaflet.num_trees
    leaflet = add_chordae(leaflet, tree_idx, chordae_copy(tree_idx).C); 
end 

valve.leaflets(1) = leaflet; 


if debug 
    load('mitral_no_partition_256_final_data.mat', 'valve_with_reference'); 
    
    tol = eps; 
    
    % [common, d1, d2] = comp_struct(valve,valve_with_reference,2,0,tol)
    
    
    [common, d1, d2] = comp_struct(valve.leaflets(1),valve_with_reference.leaflets(1),2,0,tol)
    
end 





