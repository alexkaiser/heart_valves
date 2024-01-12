function valve_with_reference = aortic_dialate_annulus(valve_with_reference)


if isfield(valve_with_reference, 'dilation_dist')
    dilation_dist = valve_with_reference.dilation_dist; 
else 
    error('must provide dilation_dist if valve.dilate_graft is true'); 
end 


leaflet = valve_with_reference.leaflets(1); 

j_max  = leaflet.j_max; 
k_max  = leaflet.k_max; 
X      = leaflet.X; 

for j=1:j_max 
    for k=1:k_max 
       [th, r, z] = cart2pol(X(1,j,k), X(2,j,k), X(3,j,k));
       [x_tmp, y_tmp, z_tmp] = pol2cart(th,r + dilation_dist,z); 
       X(:,j,k) = [x_tmp; y_tmp; z_tmp]; 
    end 
end 

valve_with_reference.leaflets(1).X = X; 

valve_with_reference.skeleton.r = valve_with_reference.skeleton.r + dilation_dist; 
valve_with_reference.skeleton.r_of_z = @(z) valve_with_reference.skeleton.r_of_z(z) + dilation_dist; 

debug = true; 
if debug 
%     figure; 
% 
%     x_component = squeeze(X(1,:,:)); 
%     y_component = squeeze(X(2,:,:)); 
%     z_component = squeeze(X(3,:,:)); 
% 
%     width = 1.0; 
%     surf(x_component, y_component, z_component, 'LineWidth',width);
%     axis equal 
%     axis auto 

    valve_plot(valve_with_reference); 
    
end 


