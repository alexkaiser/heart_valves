function Rot = rotation_matrix(angle)
% returns counter clockwise rotation matrix by angle 

    Rot = [cos(angle) sin(angle); -sin(angle) cos(angle)]; 
end 

