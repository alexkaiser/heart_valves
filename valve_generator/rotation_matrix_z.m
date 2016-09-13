function R = rotation_matrix_z(theta)
    % 
    % Rotation matrix counter clockwise by theta 
    % around z axis 

    R = [cos(theta) -sin(theta) 0; 
         sin(theta)  cos(theta) 0; 
         0           0          1]; 
end 
