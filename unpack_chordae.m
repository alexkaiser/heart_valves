function [C_left, C_right, left_papillary, right_papillary] = unpack_chordae(chordae)

C_left            = chordae.C_left; 
C_right           = chordae.C_right;  
left_papillary    = chordae.left_papillary; 
right_papillary   = chordae.right_papillary; 

[M total_len] = size(C_left); 

if M ~= 3
    error('general chordae must be 3d'); 
end 

% N = (total_len + 1) / 2; 

