function [] = tree_plot(leaflet, tree_idx, fig)
% 
% Plots chordae tendineae on figure 'fig'
% 
% leaflet     Big data strucure 
% fig         Current figure handle 

set(0, 'currentfigure', fig); 
hold on 

% unpack some data 
X              = leaflet.X; 
root           = leaflet.chordae(tree_idx).root; 
free_edge_idx  = leaflet.chordae(tree_idx).free_edge_idx; 
C              = leaflet.chordae(tree_idx).C; 


if isfield(leaflet, 'reflect_x') && leaflet.reflect_x
    X(1,:,:) = -X(1,:,:);  
    C(1,:)   = -C(1,:); 
    root(1)  = -root(1); 
end 


[n_leaves, m] = size(free_edge_idx);  
n_tree = log2(n_leaves);

% there are max_internal points on the tree 
% leaves are not included as they are part of the leaflet 
[m max_internal] = size(C); 

if m ~= 3
    error('chordae must be three dimensional')
end 

% this parameter is the same N as in the leaflet 
% it is the number of (not included) leaves in the tree
if n_leaves ~= (max_internal + 1)
    error('inconsistent dimensions of chordae and free edge indices')
end 

% sanity check in building a balanced tree 
if abs(n_tree - round(n_tree)) > eps 
    error('must use a power of two'); 
end 

width = 1; 


% left side 
x = [root(1); C(1,1)]; 
y = [root(2); C(2,1)]; 
z = [root(3); C(3,1)];
plot3(x,y,z,'k' ,'LineWidth',width); 


for i = 1:max_internal
    
    left_child_idx = 2*i; 
    
    if left_child_idx <= max_internal
        left_child = C(:,left_child_idx); 
    else 
        j = left_child_idx - max_internal; 
        left_child = X(:, free_edge_idx(j,1), free_edge_idx(j,2));
    end 
    
    x = [C(1,i); left_child(1)]; 
    y = [C(2,i); left_child(2)]; 
    z = [C(3,i); left_child(3)]; 
    plot3(x,y,z,'k' ,'LineWidth',width); 
    
   
    right_child_idx = 2*i + 1; 
    
    if right_child_idx <= max_internal
        right_child = C(:,right_child_idx); 
    else 
        j = right_child_idx - max_internal; 
        right_child = X(:, free_edge_idx(j,1), free_edge_idx(j,2));
    end 
    
    x = [C(1,i); right_child(1)]; 
    y = [C(2,i); right_child(2)]; 
    z = [C(3,i); right_child(3)];
    plot3(x,y,z,'k' ,'LineWidth',width); 
    
end 


axis auto 
axis equal 



















