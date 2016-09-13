function [C_left, C_right, left_papillary, right_papillary, Ref_l, Ref_r, k_l, k_r, k_0, k_multiplier] = unpack_chordae(chordae)

C_left            = chordae.C_left; 
C_right           = chordae.C_right;  
left_papillary    = chordae.left_papillary; 
right_papillary   = chordae.right_papillary;
Ref_l             = chordae.Ref_l; 
Ref_r             = chordae.Ref_r;  
k_l               = chordae.k_l; 
k_r               = chordae.k_r; 
k_0               = chordae.k_0; 
k_multiplier      = chordae.k_multiplier; 

[M total_len] = size(C_left); 

if M ~= 3
    error('general chordae must be 3d'); 
end 

% N = (total_len + 1) / 2; 

