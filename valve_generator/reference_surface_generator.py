'''
Generates initial and reference surfaces.

Alexander Kaiser 
9/2016
'''

import numpy as np 
import scipy
import scipy.integrate 


class cone_filter:
    '''
    Parameters for generating a coffee cone filter surface.
    
    a   Spacing of base points
    r   Radius at ring
    h   Height of cone
    N   Number of points on a leaflet
    '''
    
    def __init__(self, a, r, h, N, min_angle, max_angle):

        self.a = a
        self.r = r
        self.h = h
        self.N = N
        self.min_angle = min_angle
        self.max_angle = max_angle

    def compute_intersection(self, X, j, k):
        '''
        
        '''
        assert False 


    def apply(self, X_flat, j, k, back_sheet=False):
        '''
        Maps a position in the preimage to three dimensional image 
        '''
        assert False 
      
      
    def apply_inverse(self, X, back_sheet=False):
        ''' Maps position X to its preimage in 2d '''
        
        assert X[2] >= 0.0, 'z coordinate required to be nonnegative'
        assert not back_sheet, 'Explicit back sheet not yet implemented'
        
        tol = 1.0e-12;
        
        # maximum angle in both coordinate systems 
        theta_0 = np.arccos(-self.a / np.sqrt(self.a**2 + self.r**2 + self.h**2));

        integrand = lambda phi: (-np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2) 
                              / (1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2))

        theta_min = theta_0 + scipy.integrate.quad(integrand, 0, np.pi, epsabs=tol, epsrel=tol)
        
        psi = scipy.arcsin(self.r / np.sqrt(self.h**2 + self.r**2))
        rotated = rotation_matrix_Y(psi) * X 

        # Are we on the triangle? 
        # If so the map is a rotation 
        if abs(rotated(1)) < tol: 
            xi  = rotated(2)
            eta = rotated(3)

            if eta > 0.0: 
                theta_right = np.atan2(eta, xi - self.a) 
                assert theta_right < theta_0, 'we are supposed to be in the triangle, but are on the right cone' 
                    
                theta_left = atan2(eta, -xi - a);
                assert theta_left < theta_0, 'we are supposed to be in the triangle, but are on the left cone' 

            return np.array([xi,eta])

        # Are we on the right side cone?
        # Positive x,y  
        elif (0.0 <= X[0]) and (0.0 <= X[1]): 
            return self.apply_right_cone_inverse(X)

        # If x is negative, y positive, we have wrapped to the back sheet 
        # Reflect x, and apply front sheet map  
        elif (X[0] < 0.0) and (0.0 <= X[1]):             
            
            X[0] = -X[0]            
            val = self.apply_right_cone_inverse(X)

            # map the right front preimage to the right back preimage 
            val[0] =  out[0] - a; 
            val    =  rotation_matrix_2d(theta_min - pi/2) * val; 
            val[0] = -val[0]; 
            val    =  rotation_matrix_2d(-(theta_min - pi/2)) * val; 
            val[0] =  val[0] + a;

            return val 
            
        # left cone
        elif (0.0 <= X[0]) and (X[1] < 0.0):  
            
            # reflect and apply right cone inverse 
            X[1] = -X[1]
            val  = self.apply_right_cone_inverse(X)
            
            # back to the left 
            val[1] = -val[1]
            return val 
            
        elif (X[0] < 0.0) and (X[1] < 0.0):             
        
            # reflect back to the left front 
            X[0] = -X[0];
            
            # now use the left front inverse
            X[1] = -X[1]; 
            val = cone_filter_inv_right_cone(X, filter_params); 
            val[0] = -val[0]; 
            
            # map the left front preimage to the left back preimage 
            val[0] = val[0] + a; 
            val = rotation_matrix(-(theta_min - pi/2)) * val; 
            val[0] = -val[0]; 
            val = rotation_matrix(theta_min - pi/2) * val; 
            val[0] = val[0] - a; 
        
            return val 
        
        else: 
            assert False, 'No valid parameter range found, out of range or bugs' 
            
            
            

    def apply_right_cone_inverse(self, X): 
        assert False 

          



def rotation_matrix_Y(theta):
    ''' Returns 3x3 matrix that rotates by theta counter clockwise around z axis'''
    
    R = np.matrix(np.identity(3))
    
    c = np.cos(theta)
    s = np.sin(theta)    
    
    R[0,0] =  c 
    R[0,2] = -s 
    R[2,0] =  s 
    R[2,2] =  c 
    
    return R 
    
    
def rotation_matrix_2d(theta):
    ''' Returns 3x3 matrix that rotates by theta counter clockwise around z axis'''
    
    R = np.matrix(np.identity(2))
    
    c = np.cos(theta)
    s = np.sin(theta)    
    
    R[0,0] =  c 
    R[0,1] = -s 
    R[1,0] =  s 
    R[1,1] =  c 
    
    return R 
    
    
    
