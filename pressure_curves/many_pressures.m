% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 




cycle_length = 0.8; 

% quadrature spacing 
dt = 1e-5; 

plots = true; 

base_name = 'fourier_coeffs';


basic = false; 
if basic 

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120; 
    0.75, 130; 
    cycle_length, 8]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 20; 
    0.53, 5; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = ''; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 

low = true; 
if low 

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120*.5; 
    0.75, 130*.5; 
    cycle_length, 8]; 

    depressed_low_pressure = false; 
    if depressed_low_pressure
        points_one_cycle_atrium = [0.0, 12.5; 
        0.06, 4; 
        0.40, 7; 
        0.47, 20; 
        0.53, 5; 
        0.58, 7; 
        0.7,  8; 
        cycle_length, 12.5];
    else 
        points_one_cycle_atrium = [0.0, 24.555; 
        0.06, 4; 
        0.40, 7; 
        0.47, 20; 
        0.53, 5; 
        0.58, 7; 
        0.7,  10; 
        cycle_length, 24.555]; 
    end 
    
    suffix = '_low'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 

end 


high = false; 
if high 

    points_one_cycle_ventricle = [0.0,   -4; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 14; 
    0.58, 120 * 2; 
    0.75, 130 * 2; 
    cycle_length, -4]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 20; 
    0.53, 5; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = '_high'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 


no_atrial_systole = false; 
if no_atrial_systole

    points_one_cycle_ventricle = [0.0,   0; 
    0.02, -4; 
    0.06, 2; 
    0.40, 6; 
    0.53, 6; 
    0.58, 120; 
    0.75, 130; 
    cycle_length, 8]; 

    points_one_cycle_atrium = [0.0, 24.555; 
    0.06, 4; 
    0.40, 7; 
    0.47, 7; 
    0.53, 7; 
    0.58, 7; 
    0.7,  10; 
    cycle_length, 24.555]; 

    suffix = 'no_atrial_systole'; 

    ventricular_pressure_yellin(cycle_length, dt, points_one_cycle_ventricle, points_one_cycle_atrium, base_name, suffix); 
end 
