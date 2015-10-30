function C = cross_matrix(x)
%
% Returns the 3x3 matrix which has the action of applying x cross 
% 
    C = [   0   x(3) -x(2); 
         -x(3)    0   x(1); 
          x(2) -x(1)    0];
      
end 