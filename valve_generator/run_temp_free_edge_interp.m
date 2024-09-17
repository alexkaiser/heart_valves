function [free_edge_len, free_edge_interp_points] = run_temp_free_edge_interp(leaflet, extra_stretch_radial, y_max_from_center, power)

% runs the interpolation via input curve 
% and returns the resulting length 


k_max  = leaflet.k_max; 
N_each = leaflet.N_each; 

R_u    = leaflet.R_u; 
R_v    = leaflet.R_v; 

if isfield(leaflet, 'N_leaflets')
    N_leaflets = leaflet.N_leaflets; 
else 
    error("must provide N_leaflets")
end 

if N_leaflets ~= 2
    error("N_leaflets must be equal to 2")
end 

if ~exist('extra_stretch_radial', 'var')
    extra_stretch_radial = 1.0; 
end

if ~exist('power', 'var')
    power = 4.0; 
end

X = leaflet.X; 

% start with initial free edge 
% commissure points stay fixed 
free_edge_interp_points = X(:,:,k_max); 


for j=2:N_each

    k = k_max; 

    ring_point = X(:,j,1); 

    th = atan2(ring_point(2), ring_point(1));  

    % make a litle 
    strained_len_total = extra_stretch_radial * sum(R_v(j, :)); 

    % cm apart at middle 
%         y_free_edge_end = y_max_from_center * sign(sin(th)) * sin(th)^2; 
    % power = 4; 
    y_free_edge_end = y_max_from_center * sign(sin(th)) * abs(sin(th)^power);
    % y_free_edge_end = y_max_from_center * ring_point(2);
    % y_free_edge_end = 0; 
    % this would put the two free edges exactly coinciding 

    % if using exact x 
    % then (y_diff^2 + height^2) = strained_len_total^2 
    % so height is given as 
    interp_height = sqrt(strained_len_total^2 - (ring_point(2) - y_free_edge_end)^2) ; 

    % comm_interp_point = [ring_point(1) ; (r/2) * sin(th); comm_prev(3)];                 
    free_edge_interp_points(:,j) = [ring_point(1) ; y_free_edge_end; interp_height + ring_point(3)]; 

    % total radial rest length of this radial fiber 
    total_rest_length = sum(R_v(j, 1:(k-1))); 

    tangent = (free_edge_interp_points(:,j) - ring_point); 
    tangent = tangent / norm(tangent); 

    % based on the rest length 
    X(:,j,k) = total_rest_length * tangent * extra_stretch_radial + ring_point; 

end 

free_edge_len = get_circ_edge_lengths(leaflet, N_each, k_max, X, R_u); 
