function [free_edge_idx_left free_edge_idx_right] = get_free_edge_ranges(leaflet)
% 
% Returns two 2d arrays of indices 
% 
% Input: 
%    leaflet    Parameter struct 
% 
% Output: 
%    free_edge_idx_left    2d array of left  chordae indices 
%    free_edge_idx_right   2d array of right chordae indices  
%

if leaflet.radial_and_circumferential
    error('radial and circumferential not implemented')
end 

N = leaflet.N; 

free_edge_idx_left  = [ones(N,1), (1:N)'];
free_edge_idx_right = [(1:N)', ones(N,1)];
