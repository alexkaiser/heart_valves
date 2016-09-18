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
R       = leaflet.R; 
chordae = leaflet.chordae; 

if left_side 
    C     = chordae.C_left; 
    R_ch  = chordae.Ref_l; 
    pap   = chordae.left_papillary; 
    k_spr = chordae.k_l; 
    k_0   = chordae.k_0; 
else 
    C     = chordae.C_right; 
    R_ch  = chordae.Ref_r; 
    pap   = chordae.right_papillary;
    k_spr = chordae.k_r; 
    k_0   = chordae.k_0; 
end

[m max_internal] = size(C); 

% parent direction neighbor may be the papillary muscle
% this occurs precisely when requesting the zero index
if nbr_idx == 0 
    nbr   = pap; 
    R_nbr = pap; 

% if neighbors are out of the chordae region
% then they may be on the leaflet 
elseif nbr_idx > max_internal
    
    free_edge_idx = nbr_idx - max_internal; 
    if left_side 
        j = leaflet.free_edge_idx_left (free_edge_idx,1); 
        k = leaflet.free_edge_idx_left (free_edge_idx,2); 
    else 
        j = leaflet.free_edge_idx_right(free_edge_idx,1); 
        k = leaflet.free_edge_idx_right(free_edge_idx,2); 
    end 
    
    nbr   = X(:,j,k); 
    R_nbr = R(:,j,k); 

% the neighbor is within the tree of chordae 
else 
    nbr   = C(:,nbr_idx); 
    R_nbr = R_ch(:,nbr_idx); 
end 

% spring constants
if nbr_idx < i 
    % nbr_idx is only less if nbr is the parent 
    % parent wise owned by this index 
    k_val = k_spr(i);       

% connections from chordae to leaflet have value k_0
elseif nbr_idx > max_internal 

    k_val = k_0;                 

else 

    % child's parent-direction spring is at child's index 
    k_val = k_spr(nbr_idx);     

end 




