function eta = compute_elastic_modulus(valve)
% 
% Computes the approximate elastic modulus eta for the given set of paramters 
% 


    
N                            = valve.N; 
L                            = valve.L; 
pressure_tension_ratio       = valve.pressure_tension_ratio; 
p_physical                   = valve.p_physical; 
refinement                   = valve.refinement; 


% approximate length element 
ds = 2*L / N;

% physical thickness 
thickness = 0.05; 

dA = ds * thickness; 

MMHG_TO_CGS = 1333.22368; 
p_cgs = p_physical * MMHG_TO_CGS; 

% value of the relative spring constant is determined by the ratio 
k_rel = p_cgs / pressure_tension_ratio; 

% relative spring constants drop when the mesh is refined 
k_rel = k_rel / refinement; 

eta = k_rel / dA; 







