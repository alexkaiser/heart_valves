function [nbr R_nbr k_val j k] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx)
% 
% Given a current index and nieghbor index 
% Returns the coordinates, reference coordinates and spring constant
% Takes into account all boundary conditions 
% 
% Input
%     leaflet     Current main data structure 
%     i           Index in the chordae tree 
%     nbr_idx     Index of the neighbor in chordae tree 
%     left_side   True if on the left tree 
% 
% Output 
%     nbr         Coordinates of the neighbor
%     R_nbr       Reference coordinate of the neighbor
%     k_val       Spring constant of the connector 
%     j,k         If the neighbor is on the leaflet, these are its coordinates 
%                 Empty if the neighbor is not on the leaflet 


% default empty values 
j     = []; 
k     = []; 
R_nbr = []; 

X       = leaflet.X; 
chordae = leaflet.chordae; 

if isfield(chordae(tree_idx), 'R_free_edge')  && isfield(chordae(tree_idx), 'k_free_edge')
    
    free_edge_constants_set = true; 
    
    R_free_edge = chordae(tree_idx).R_free_edge;
    k_free_edge = chordae(tree_idx).k_free_edge;
    
else
    
    free_edge_constants_set = false;

end 

C             = chordae(tree_idx).C; 
root          = chordae(tree_idx).root; 
k_vals        = chordae(tree_idx).k_vals; 
k_0           = chordae(tree_idx).k_0; 
free_edge_idx = chordae(tree_idx).free_edge_idx; 

if isfield(chordae(tree_idx), 'R_ch')
    R_ch  = chordae(tree_idx).R_ch; 
else 
    R_ch  = []; 
end 


[m max_internal] = size(C); 

% if neighbors are out of the chordae region
% then they may be on the leaflet 
if nbr_idx > max_internal
    
    idx = nbr_idx - max_internal; 

    j = free_edge_idx(idx,1); 
    k = free_edge_idx(idx,2); 
    
    nbr   = X(:,j,k); 
        
    if free_edge_constants_set
        % fetch from free edge arrays if available 
        R_nbr = R_free_edge(idx); 
        k_val = k_free_edge(idx); 
    else 
        % free edge springs are all k_0 if not  
        R_nbr = []; 
        k_val = k_0;
    end 
    
% the neighbor is within the tree of chordae 
else 
    
    % parent direction neighbor may be the papillary muscle
    % this occurs precisely when requesting the zero index
    if nbr_idx == 0 
        
        nbr   = root; 
        
        % papillary muscle is always parent, so current location owns constant
        if ~isempty(R_ch)
            R_nbr = R_ch(i); 
        end
        
        k_val = k_vals(i);
        
    else

        nbr   = C(:,nbr_idx); 
        
        % spring constants
        if nbr_idx < i 
            % nbr_idx is only less if nbr is the parent 
            % parent wise owned by this index 
            k_val = k_vals(i);                      
            
            if ~isempty(R_ch)
                R_nbr = R_ch(i); 
            end
            
        else 

            % child's parent-direction spring is at child's index 
            k_val = k_vals(nbr_idx);
            
            if ~isempty(R_ch)
                R_nbr = R_ch(nbr_idx); 
            end 
        end 
    
    end 
end 

