function fig = dissection_plot_rest_height(valve, fig)

if ~exist('fig', 'var')
    fig = figure; 
end 

if isfield(valve, 'name') && strcmp(valve.name, 'aortic') 
    warning('dissection_plot_rest_height not implemented for aortic'); 
    return; 
end 
    

valve_with_reference = valve; 
valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 

valve_with_reference.leaflets(1) = set_rest_lengths_and_constants(valve.leaflets(1), valve); 

leaflet = valve_with_reference.leaflets(1); 

X_current              = leaflet.X; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
R_v                    = leaflet.R_v; 

circ_cumulative_sum = zeros(1,j_max); 

circ_cumulative_sum(1) = 0.0; 

for j=2:j_max
    circ_cumulative_sum(j) = circ_cumulative_sum(j-1) + norm(X_current(:,j-1,k_max) - X_current(:,j,k_max)); 
end 

positions_x = nan(j_max,k_max); 
positions_y = nan(j_max,k_max); 

for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) || is_bc(j,k)
            positions_x(j,k) = circ_cumulative_sum(j); 
            
            positions_y(j,k) = 0;

            if k<k_max 
                for k_tmp = k:(k_max-1)
                    
                    j_nbr_tmp = j; 
                    k_nbr_tmp = k_tmp + 1; 
                    
                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k_tmp, j_nbr_tmp, k_nbr_tmp); 

                    if ~valid 
                        error('trying to compute lengths with an invalid rest length')
                    end 
                    positions_y(j,k) = positions_y(j,k) - R_v(j_spr,k_spr); 

                end 
                
            end 
            
        end 
    end 
end 

plot(positions_x, positions_y, 'k.'); 
axis equal
xlabel('circumferential (cm)')
ylabel('radial (cm)')
title('radial height only')


j_center_anterior  = leaflet.j_range_anterior(floor(end/2)); 
j_center_posterior = leaflet.j_range_posterior(floor(end/2)); 

k_center_anterior = find(~isnan(positions_y(j_center_anterior,:)), 1); 
height_anterior = abs(positions_y(j_center_anterior, k_center_anterior)); 

k_center_posterior = find(~isnan(positions_y(j_center_posterior,:)), 1); 
height_posterior = abs(positions_y(j_center_posterior, k_center_posterior)); 

fprintf('Rest length height summary:\n'); 
fprintf('Anterior  height = %.3f cm\n', height_anterior)
fprintf('Posterior height = %.3f cm\n', height_posterior)








