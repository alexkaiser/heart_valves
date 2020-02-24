function fig = dissection_plot_rest_height_aortic(valve, fig)

if ~exist('fig', 'var')
    fig = figure; 
end 

if ~isfield(valve, 'name') || ~strcmp(valve.name, 'aortic') 
    error('dissection_plot_rest_height_aortic called with non aortic type'); 
    return; 
end 
    
valve_with_reference = valve; 
valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 

valve_with_reference.leaflets(1) = set_rest_lengths_and_constants_aortic(valve.leaflets(1), valve); 

leaflet = valve_with_reference.leaflets(1); 

radius = valve.skeleton.r; 

X_current              = leaflet.X; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
R_v                    = leaflet.R_v; 
R_u                    = leaflet.R_u; 
N_each                 = leaflet.N_each; 

circ_cumulative_sum = linspace(0, 2*pi*radius, j_max+1); 
circ_cumulative_sum = circ_cumulative_sum(1:j_max);  

positions_x = nan(j_max,k_max); 
positions_y = nan(j_max,k_max); 

center_leaflet_idx = N_each/2; 

for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) || is_bc(j,k)
            positions_x(j,k) = circ_cumulative_sum(j); 
            
            positions_y(j,k) = X_current(3,j,1); % z component at bottom 

            if j == center_leaflet_idx 
                'stop'; 
            end 
            
            if k > 1 
                for k_tmp = 1:(k-1)
                    
                    j_nbr_tmp = j; 
                    k_nbr_tmp = k_tmp + 1; 
                    
                    [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k_tmp, j_nbr_tmp, k_nbr_tmp); 

                    if ~valid 
                        error('trying to compute lengths with an invalid rest length')
                    end 
                    positions_y(j,k) = positions_y(j,k) + R_v(j_spr,k_spr); 

                end 
                
            end 
            
        end 
    end 
end 

subplot(2,1,1)
plot(positions_x, positions_y, 'k.'); 
axis equal
xlabel('circumferential (cm)')
ylabel('radial (cm)')
title('aortic, radial height only')



length_one_leaflet_free_edge = zeros(k_max,1); 
for j=1:N_each
    for k=1:k_max
        if is_internal(j,k) || is_bc(j,k)
            j_nbr_tmp = j+1; 
            k_nbr_tmp = k; 
            [valid j_nbr k_nbr j_spr k_spr target_spring] = get_indices(leaflet, j, k_tmp, j_nbr_tmp, k_nbr_tmp); 

            if ~valid 
                error('trying to compute lengths with an invalid rest length')
            end

            length_one_leaflet_free_edge(k) = length_one_leaflet_free_edge(k) + R_u(j_spr,k_spr); 
        end 
    end 
end 

subplot(2,1,2)
j = center_leaflet_idx; 
for k=1:k_max
    y = [positions_y(j,k) positions_y(j,k)]; 
    x = length_one_leaflet_free_edge(k)/2 * [-1, 1]; 
    plot(x,y,'k'); 
    hold on 
end 
    
% r = d/2 (clearly)
% h = 0.71d = 1.42r, height to commissure 
% c = .17d = .34r, coaptation height
% ?f = .62 = 1.24r, free edge(is this entire free edge length?). this must be half length. 
% lc = .7d = 1.4r. the inner length of the dotted line on the left,  leaflet height in 2d. 

fprintf('Rest length height summary:\n'); 
fprintf("and targets based on Swanson and Clark 1973\n")
fprintf("Radius                  = %f\n", radius); 
fprintf("One third circumference = %f\n", 2*pi*radius/3); 
fprintf("Free edge target        = %f\n", 2*1.24*radius); 
fprintf("Circ free edge length   = %f\n", length_one_leaflet_free_edge(k_max)); 
fprintf("Leaflet height target   = %f\n", 1.4*radius); 
fprintf("Radial height           = %f\n", radial_leaflet_height); 

% 
% j_center_anterior  = leaflet.j_range_anterior(floor(end/2)); 
% j_center_posterior = leaflet.j_range_posterior(floor(end/2)); 
% 
% k_center_anterior = find(~isnan(positions_y(j_center_anterior,:)), 1); 
% height_anterior = abs(positions_y(j_center_anterior, k_center_anterior)); 
% 
% k_center_posterior = find(~isnan(positions_y(j_center_posterior,:)), 1); 
% height_posterior = abs(positions_y(j_center_posterior, k_center_posterior)); 
% 
% fprintf('Rest length height summary:\n'); 
% fprintf('Anterior  height = %.3f cm\n', height_anterior)
% fprintf('Posterior height = %.3f cm\n', height_posterior)
% 
% 






