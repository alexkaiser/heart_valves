function [] = tree_plot(leaflet, fig)
% 
% Plots chordae tendineae on figure 'fig'
% 
% params      Big data strucure 
% fig         Current figure handle 

set(0, 'currentfigure', fig); 
hold on 

% unpack some data 
X                   = leaflet.X; 
left_papillary      = leaflet.left_papillary;
right_papillary     = leaflet.right_papillary;
free_edge_idx_left  = leaflet.free_edge_idx_left; 
free_edge_idx_right = leaflet.free_edge_idx_right; 
C_left              = leaflet.chordae.C_left; 
C_right             = leaflet.chordae.C_right; 

if leaflet.reflect_x
    X(1,:,:)           = -X(1,:,:);  
    C_left(1,:)        = -C_left(1,:); 
    C_right(1,:)       = -C_right(1,:); 
    left_papillary(1)  = -left_papillary(1); 
    right_papillary(1) = -right_papillary(1);
end 


[n_leaves, m] = size(leaflet.free_edge_idx_left);  
n_tree = log2(n_leaves);

% there are max_internal points on the tree 
% leaves are not included as they are part of the leaflet 
[m max_internal] = size(C_left); 

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

width = 1.5; 


% left side 
x = [left_papillary(1); C_left(1,1)]; 
y = [left_papillary(2); C_left(2,1)]; 
z = [left_papillary(3); C_left(3,1)];
plot3(x,y,z,'k' ,'LineWidth',width); 

x = [right_papillary(1); C_right(1,1)]; 
y = [right_papillary(2); C_right(2,1)]; 
z = [right_papillary(3); C_right(3,1)]; 
plot3(x,y,z,'k' ,'LineWidth',width); 


for i = 1:max_internal
    
    left_child_idx = 2*i; 
    
    % left side chordae takes a k index when hitting leaflet 
    if left_child_idx <= max_internal
        left_child = C_left(:,left_child_idx); 
    else 
        k = left_child_idx - max_internal; 
        left_child = X(:, free_edge_idx_left(k,1), free_edge_idx_left(k,2));
    end 
    
    x = [C_left(1,i); left_child(1)]; 
    y = [C_left(2,i); left_child(2)]; 
    z = [C_left(3,i); left_child(3)]; 
    plot3(x,y,z,'k' ,'LineWidth',width); 
    
    % right side chordae takes a j index when hitting leaflet 
    if left_child_idx <= max_internal
        left_child = C_right(:,left_child_idx); 
    else
        j = left_child_idx - max_internal; 
        left_child = X(:, free_edge_idx_right(j,1), free_edge_idx_right(j,2));
    end 
    
    x = [C_right(1,i); left_child(1)]; 
    y = [C_right(2,i); left_child(2)]; 
    z = [C_right(3,i); left_child(3)]; 
    plot3(x,y,z,'k' ,'LineWidth',width); 
    
    
    right_child_idx = 2*i + 1; 
    
    if right_child_idx <= max_internal
        right_child = C_left(:,right_child_idx); 
    else 
        k = right_child_idx - max_internal; 
        right_child = X(:, free_edge_idx_left(k,1), free_edge_idx_left(k,2));
    end 
    
    
    x = [C_left(1,i); right_child(1)]; 
    y = [C_left(2,i); right_child(2)]; 
    z = [C_left(3,i); right_child(3)];
    plot3(x,y,z,'k' ,'LineWidth',width); 
    
    if right_child_idx <= max_internal
        right_child = C_right(:,right_child_idx); 
    else 
        j = right_child_idx - max_internal; 
        right_child = X(:, free_edge_idx_right(j,1), free_edge_idx_right(j,2));
    end 
    
    x = [C_right(1,i); right_child(1)]; 
    y = [C_right(2,i); right_child(2)]; 
    z = [C_right(3,i); right_child(3)];
    plot3(x,y,z,'k' ,'LineWidth',width); 
    % axis(axes_vec); 
    
end 


axis auto 
axis equal 



















