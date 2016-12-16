function [E, F, J] = energy_gradient_hessian(X_linearized, leaflet)
% 
% Returns energy, gradient and hessian for current leaflet data structures 
% 
% 


if ~isfield(leaflet, 'energy')
    error('Cannot run optimization without energy'); 
end

leaflet_copy = internal_points_to_2d(X_linearized, leaflet);

E = leaflet.energy(leaflet_copy); 

% Difference equations are gradient of energy 
if nargout > 1
    [F_2d F_chordae_left F_chordae_right] = leaflet.diff_eqns(leaflet_copy); 
    
    % Note that force is NEGATIVE gradient, but optimization just takes the gradient 
    F = -linearize_internal_points(leaflet, F_2d, F_chordae_left, F_chordae_right);
end

% Jacobian of gradient, Hessian of energy 
% Similarly, this is the second deriv of energy, so is NEGATIVE Jacobian of force
if nargout > 2
    J = -leaflet.jacobian(leaflet_copy); 
end


if isfield(leaflet, 'repulsive_potential') && leaflet.repulsive_potential

    E = E + energy_repulsive(leaflet_copy); 

    % Difference equations are gradient of energy 
    if nargout > 1
        [F_2d F_chordae_left F_chordae_right] = difference_equations_repulsive(leaflet_copy); 

        % Note that force is NEGATIVE gradient, but optimization just takes the gradient 
        F = F + -linearize_internal_points(leaflet, F_2d, F_chordae_left, F_chordae_right);
    end

    % Jacobian of gradient, Hessian of energy 
    % Similarly, this is the second deriv of energy, so is NEGATIVE Jacobian of force
    if nargout > 2
        J = J + -build_jacobian_replusive(leaflet_copy); 
    end
    

end 
