

N = 8; 

% Initialize structures  
valve = initialize_valve_data_structures_radial_bead_slip(N); 

'orignal data'

valve.anterior.X
valve.posterior.X
valve.anterior.chordae.C_left 
valve.anterior.chordae.C_right

% 'linear order on internal points'
X_linearized = linearize_internal_points_bead_slip(valve, valve.anterior.X, valve.posterior.X, valve.anterior.chordae.C_left, valve.anterior.chordae.C_right)


% 'after return to normal data structure'
% posterior = internal_points_to_2d(X_and_chordae_linearized, posterior);  
