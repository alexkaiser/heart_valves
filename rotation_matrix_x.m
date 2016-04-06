function R = rotation_matrix_x(theta)
    % 
    % Rotation matrix counter clockwise by theta 
    % around z axis 

    R = [1          0           0; 
         0 cos(theta) -sin(theta);  
         0 sin(theta)  cos(theta)]; 
end 
