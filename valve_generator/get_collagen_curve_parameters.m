function [collagen_curve] = get_collagen_curve_parameters()

MPa_TO_CGS = 1.0e7; 

collagen_curve.a = 4643.4; 
collagen_curve.b                    = 49.9643;   
collagen_curve.full_recruitment     = 0.145;     
collagen_curve.eta_collagen         = 32.5 * MPa_TO_CGS; 
collagen_curve.collagen_x_intercept = 0.125;             
collagen_curve.collagen_y_intercept = -collagen_curve.collagen_x_intercept * collagen_curve.eta_collagen; 


