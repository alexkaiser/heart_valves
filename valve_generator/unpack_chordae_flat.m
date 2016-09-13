function [C_left_flat, C_right_flat, left_papillary_flat, right_papillary_flat] = unpack_chordae_flat(chordae_flat)

C_left_flat            = chordae_flat.C_left_flat; 
C_right_flat           = chordae_flat.C_right_flat;  
left_papillary_flat    = chordae_flat.left_papillary_flat; 
right_papillary_flat   = chordae_flat.right_papillary_flat; 

% [M total_len] = size(C_left_flat); 
% 
% if M ~= 2
%     error('flat chordae must be 2d'); 
% end 
% 
% N = (total_len + 1) / 2; 
% 
% 

