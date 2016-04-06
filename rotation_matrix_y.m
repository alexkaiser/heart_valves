
function R = rotation_matrix_y(theta)
    % 
    % Rotation matrix counter clockwise by theta 
    % around z axis 

    R = [cos(theta)  0 -sin(theta); 
         0           1  0         ; 
         sin(theta)  0  cos(theta)]; 
end 
