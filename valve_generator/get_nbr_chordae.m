function [nbr R_nbr k_val j k] = get_nbr_chordae(leaflet, i, nbr_idx, left_side)
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
j=[]; 
k=[]; 

X       = leaflet.X; 
chordae = leaflet.chordae; 

if isfield(leaflet, 'R_free_edge_left')  && isfield(leaflet, 'k_free_edge_left') && ... 
   isfield(leaflet, 'R_free_edge_right') && isfield(leaflet, 'k_free_edge_right') 
    
    free_edge_constants_set = true; 
    
    if left_side 
        R_free_edge = leaflet.R_free_edge_left;
        k_free_edge = leaflet.k_free_edge_left;
    else 
        R_free_edge = leaflet.R_free_edge_right;
        k_free_edge = leaflet.k_free_edge_right;
    end 
else
    
    free_edge_constants_set = false;

end 



if left_side 
    C     = chordae.C_left; 
    R_ch  = chordae.Ref_l; 
    pap   = chordae.left_papillary; 
    k_spr = chordae.k_l; 
    k_0   = chordae.k_0; 
    free_edge_idx = leaflet.free_edge_idx_left; 
else 
    C     = chordae.C_right; 
    R_ch  = chordae.Ref_r; 
    pap   = chordae.right_papillary;
    k_spr = chordae.k_r; 
    k_0   = chordae.k_0; 
    free_edge_idx = leaflet.free_edge_idx_right; 
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
        
        nbr   = pap; 
        
        % papillary muscle is always parent, so current location owns constant
        R_nbr = R_ch(i); 
        k_val = k_spr(i);
        
    else

        nbr   = C(:,nbr_idx); 
        
        % spring constants
        if nbr_idx < i 
            % nbr_idx is only less if nbr is the parent 
            % parent wise owned by this index 
            k_val = k_spr(i);                      
            R_nbr = R_ch(i); 
        else 

            % child's parent-direction spring is at child's index 
            k_val = k_spr(nbr_idx);     
            R_nbr = R_ch(nbr_idx); 
        end 
    
    end 
end 

