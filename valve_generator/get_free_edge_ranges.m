function [free_edge_idx_left free_edge_idx_right chordae_idx_left chordae_idx_right] = get_free_edge_ranges(leaflet)
% 
% Returns two 2d arrays of indices 
% 
% Input: 
%    leaflet    Parameter struct 
% 
% Output: 
%    free_edge_idx_left    2d array of left chordae indices of the implicitly defined leaves. 
%                          The i-th leaf is not included in the tree, 
%                          but instead on the leaflet as 
%                          X(:,free_edge_idx_left(:,1),free_edge_idx_left(:,1))  
%                          
%    free_edge_idx_right   2d array of right chordae indices  
%
%    chordae_idx_left      If chordae_idx_left(j,k) is nonzero, then this contains 
%                          the leaf index which is connected to X(:,j,k)
% 
 

if leaflet.radial_and_circumferential
    error('radial and circumferential not implemented')
end 

N = leaflet.N; 
free_edge_idx_left  = [ones(N,1), (1:N)'];
free_edge_idx_right = [(1:N)', ones(N,1)];

chordae_idx_left  = zeros(N+1,N+1); 

j=1; 
for k=1:N
    chordae_idx_left(j,k) = k; 
end 

chordae_idx_right = zeros(N+1,N+1);

k=1; 
for j=1:N
    chordae_idx_right(j,k) = j; 
end 

