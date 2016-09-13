'''
Generates initial and reference surfaces.

Alexander Kaiser 
9/2016
'''

import numpy as np 
import scipy
import scipy.integrate 
import warnings


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


    def apply(self, X_flat, back_sheet=False):
        ''' Maps a position X_flat in the preimage to three dimensional image '''
        
        assert not back_sheet
        
        tol = 1.0e-14 
            
        eta = X_flat.item(0)
        xi  = X_flat.item(1)

        # If on the pole return it
        if (abs(eta) < tol) and (abs(xi - self.a) < tol):
            return np.array([0.0, self.a, 0.0]) 
         
        # If on the other pole return it
        if (abs(eta) < tol) and (abs(xi + self.a) < tol): 
            return np.array([0.0, self.a, 0.0]) 

        # maximum angle in both coordinate systems 
        theta_0 = np.arccos(-self.a / np.sqrt(self.a**2 + self.r**2 + self.h**2))

        # minimum angle in both coordinate systems to be on the front sheet
        integrand = lambda phi: (- np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2) 
                                 / (1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2) ) 
    
        val, err_est = scipy.integrate.quadrature(integrand, 0, np.pi/2.0, tol=tol, rtol=tol)

        if err_est > tol: 
            raise ValueError('Quadrature error greater than tolerance')
        
        theta_min = theta_0 + val 
    
        # right polar coordinate angle first 
        theta_right = np.arctan2(eta, xi - self.a)

        # Other angle 
        # use reflected coordinates, (-xi,eta)
        # and the (a,0) polar coordinate origin 
        theta_left = np.arctan2(eta, -xi - self.a)

        # are we in the triangle portion, which simply gets rotated out? 
        if (theta_0 <= theta_left) and (theta_0 <= theta_right): 

            # we are on the triangle, rotate out 
            psi = np.arcsin(self.r / np.sqrt(self.h**2 + self.r**2)) 
            return np.array([np.sin(psi) * eta, xi, np.cos(psi) * eta])   


        # If angle is less than the leftmost ray on the front 
        # then point is on the right, front sheet 
        elif (theta_min <= theta_right) and (theta_right <= theta_0): 
            return self.apply_right_cone(X_flat)
             
        # are we on the back of the right sheet 
        elif ((2*theta_min - theta_0) <= theta_right) and (theta_right <= theta_min):

            # translate, rotate and reflect to use the other coordinates 

            # reset origin 
            xi = xi - self.a 

            # rotate axis to positive eta axis
            val = rotation_matrix_2d(theta_min - np.pi/2.0) * np.matrix([[xi], [eta]])  
        
            # reflect over eta axis, negate xi 
            val[0] = -val[0] 
    
            # rotate back 
            val = rotation_matrix_2d(-(theta_min - np.pi/2.0)) * val   
            
            # unset origin
            val[0] = val[0] + self.a 
            
            # apply filter map 
            X = self.apply_right_cone(val) 
            
            # negate x for the back sheet 
            X[0] = -X[0] 
            
        # left front sheet 
        elif (theta_min <= theta_left) and (theta_left <= theta_0): 
            
            # apply the reflection in xi to get the corresponding preimage 
            val = np.array([-xi, eta])
            
            # map on rhs 
            X = self.apply_right_cone(val) 
    
            # reflect back in y in the 3d setting  
            X[1] = -X[1] 
    
        elif ((2.0*theta_min - theta_0) <= theta_left) and (theta_left <= theta_min): 
            
            # same as above but various signs are swapped 
            
            # reset origin 
            xi = xi + self.a 
    
            # rotate axis to positive eta axis
            val = rotation_matrix_2d(-(theta_min - np.pi/2)) * np.matrix([[xi], [eta]])  
             
            # reflect over eta axis, negate xi 
            val[0] = -val[0]  
    
            # rotate back 
            val = rotation_matrix_2d(theta_min - np.pi/2) * val   
            
            # unset origin
            val[0] -= self.a  
            
            # apply filter map using the LEFT map 
            val[0] = -val[0]
            X = self.apply_right_cone(val) 
            
            # negate in y for the left side 
            X[1] = -X[1] 
            
            # negate x for the back sheet 
            X[0] = -X[0]  
            
        else:
            raise ValueError('Should not have gotten here, coordinates not in any acceptable regions') 
        

        return X 
         
         
    def apply_right_cone(self, X_flat):
        ''' 
        Apply the right side cone
        
        Arguments must be in range 
        '''
        
        tol = 1.0e-14
        
        eta = X_flat.item(0)
        xi  = X_flat.item(1)
        
        # If on the pole return it
        if (abs(eta) < tol) and (abs(xi - self.a) < tol):  
            return np.array([0, self.a, 0])  
 
        # angle first 
        theta = np.arctan2(eta, xi - self.a)
        
        dphi_dtheta = lambda phi, th:  (-(1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2) 
                                    / np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2)) 
        
        theta_0 = np.arccos(-self.a / np.sqrt(self.a**2 + self.r**2 + self.h**2)) 
        
        if (theta - theta_0) > tol: 
            raise ValueError('requesting theta greater than maximum allowed')
 
        integrand = lambda phi: ( -np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2) 
                                  / (1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2))
        
        
        val, err_est = scipy.integrate.quadrature(integrand, 0, np.pi/2.0, tol=tol, rtol=tol)

        if err_est > tol: 
            raise ValueError('Quadrature error greater than tolerance')
        
        theta_min = theta_0 + val  
        
        if (theta - theta_min) < -tol:
            raise ValueError('requesting theta less than minimum allowed') 
         
        theta_span = [theta_0, theta] 
        phi_0 = 0  
        
        if abs(theta - theta_0) > tol: 
        
            # solve ODE for phi
            ode_solve_tol = 1e-10
            
            phi = scipy.integrate.odeint(dphi_dtheta, phi_0, theta_span, rtol=ode_solve_tol, atol=ode_solve_tol) 
            
            # always first element 
            assert len(phi) == 2
            phi = phi[1]
             
        else: 
            # theta is the initial condition, so we know phi here without the ODE 
            phi = phi_0 
 
        
        R = np.sqrt(self.r**2 - 2*self.a*self.r*np.sin(phi) + self.a**2 + self.h**2)  
        
        # catch for vertical points where first z formula does not work
        if abs(theta - np.pi/2) < tol: 
            print 'hit vertical catch' 
            z = self.h * eta / (R * np.sin(theta)) 
        elif (abs(theta - np.pi) < tol) or ((abs(theta) < tol)):
            raise ValueError('cannot take horizontal here') 
        else:     
            z1 = self.h * (xi - self.a) / (R * np.cos(theta)) 
            z2 = self.h * eta / (R * np.sin(theta)) 
        
            if abs(z1 - z2) > 1.0e3*tol:
                print 'two formulas for z disagree, z1 = ', z1, ', z2 = ', z2 
            else:
                z = z1  
        
        x = z*(self.r/self.h) * np.cos(phi) 
        y = self.a*(1 - z/self.h) + z*(self.r/self.h) * np.sin(phi) 
        
        return np.array([x, y, z])
 
      
    def apply_inverse(self, X, back_sheet=False):
        ''' 
        Maps position X to its preimage in 2d 
        
        
        '''
        
        if X[2] < 0.0:
            raise ValueError('z coordinate required to be nonnegative')
        
        assert not back_sheet, 'Explicit back sheet not yet implemented'
        
        tol = 1.0e-12
        
        # maximum angle in both coordinate systems 
        theta_0 = np.arccos(-self.a / np.sqrt(self.a**2 + self.r**2 + self.h**2))

        integrand = lambda phi: (-np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2) 
                              / (1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2))

        val, err_est = scipy.integrate.quadrature(integrand, 0, np.pi, tol=tol, rtol=tol)

        if err_est > tol: 
            raise ValueError('Quadrature error greater than tolerance')

        theta_min = theta_0 + val
        
        psi = scipy.arcsin(self.r / np.sqrt(self.h**2 + self.r**2))
        rotated = rotation_matrix_Y(psi) * X 

        # Are we on the triangle? 
        # If so the map is a rotation 
        if abs(rotated[0]) < tol: 
            xi  = rotated.item(1)
            eta = rotated.item(2)

            if eta > 0.0: 
                theta_right = np.arctan2(eta, xi - self.a) 
                if theta_right > theta_0:
                    raise ValueError('we are supposed to be in the triangle, but are on the right cone')
                    
                theta_left = np.arctan2(eta, -xi - self.a)
                if theta_left > theta_0:
                    raise ValueError('we are supposed to be in the triangle, but are on the left cone')

            return np.array([xi,eta])

        # Are we on the right side cone?
        # Positive x,y  
        elif (0.0 <= X[0]) and (0.0 <= X[1]): 
            val = self.apply_right_cone_inverse(X)

        # If x is negative, y positive, we have wrapped to the back sheet 
        # Reflect x, and apply front sheet map  
        elif (X[0] < 0.0) and (0.0 <= X[1]):             
            
            X[0] = -X[0]            
            val = self.apply_right_cone_inverse(X)

            # map the right front preimage to the right back preimage 
            val[0] =  val[0] - self.a 
            val    =  rotation_matrix_2d(theta_min - np.pi/2) * val 
            val[0] = -val[0] 
            val    =  rotation_matrix_2d(-(theta_min - np.pi/2)) * val 
            val[0] =  val[0] + self.a 
            
        # left front side of cone
        elif (0.0 <= X[0]) and (X[1] < 0.0):  
            
            # reflect and apply right cone inverse 
            X[1] = -X[1]
            val  = self.apply_right_cone_inverse(X)
            
            # back to the left 
            val[1] = -val[1]
 
        # left back side of cone            
        elif (X[0] < 0.0) and (X[1] < 0.0):             
        
            # reflect back to the left front 
            X[0] = -X[0]
            
            # now use the left front inverse
            X[1]   = -X[1] 
            val    = self.apply_right_cone_inverse(X) 
            val[0] = -val[0] 
            
            # map the left front preimage to the left back preimage 
            val[0] = val[0] + self.a 
            val    = rotation_matrix_2d(-(theta_min - np.pi/2)) * val 
            val[0] = -val[0] 
            val    = rotation_matrix_2d(theta_min - np.pi/2) * val 
            val[0] = val[0] - self.a 
               
        else: 
            raise ValueError('No valid parameter range found, out of range or bugs') 
            
        return np.array([val.item(0), val.item(1)])
          

    def apply_right_cone_inverse(self, X): 
        ''' Apply the inverse of the right cone map '''
        
        tol = 1.0e-14
        
        x = X.item(0)
        y = X.item(1)
        z = X.item(2)
        
        if z <= 0.0:
            raise ValueError('Cannot evaluate cone for negative z') 

        phi_1 = np.arccos(  x           * self.h/(z*self.r))
        phi_2 = np.arcsin( (y - self.a) * self.h/(z*self.r) + self.a/self.r) 
        if abs(phi_1 - phi_2) > tol:  
            print 'phi_1 = ', phi_1, 'phi_2 = ', phi_2, ' diff = ', abs(phi_1 - phi_2)
            warnings.warn('two formulas for phi disagree')
            phi = phi_1 
        else:
            phi = phi_1

        if (phi < 0) or (phi > np.pi/2):
            raise ValueError('phi out of range')

        R = np.sqrt(self.r**2 - 2*self.a*self.r*np.sin(phi) + self.a**2 + self.h**2) 

        theta_0 = np.arccos(-self.a / np.sqrt(self.a**2 + self.r**2 + self.h**2)) 

        integrand = lambda phi: (- np.sqrt( (1 - (self.a/self.r) * np.sin(phi))**2 + (self.h/self.r)**2) 
                                    / (1 - (2*self.a/self.r) * np.sin(phi) + (self.a/self.r)**2 + (self.h/self.r)**2))  

        val, err_est = scipy.integrate.quadrature(integrand, 0, phi, tol=tol, rtol=tol)

        if err_est > tol: 
            raise ValueError('Quadrature error greater than tolerance')

        theta = theta_0 + val

        xi  = self.a + (R*z/self.h) * np.cos(theta) 
        eta =          (R*z/self.h) * np.sin(theta) 

        return np.matrix(np.array([[xi], [eta]])) 



def rotation_matrix_Y(theta):
    ''' Returns 3x3 matrix that rotates by theta counter clockwise around y axis'''
    
    R = np.matrix(np.identity(3))
    
    c = np.cos(theta)
    s = np.sin(theta)    
    
    R[0,0] =  c 
    R[0,2] = -s 
    R[2,0] =  s 
    R[2,2] =  c 
    
    return R 
    
    
def rotation_matrix_2d(theta):
    ''' Returns 2x2 matrix that rotates by theta counter clockwise'''
    
    R = np.matrix(np.identity(2))
    
    c = np.cos(theta)
    s = np.sin(theta)    
    
    R[0,0] =  c 
    R[0,1] = -s 
    R[1,0] =  s 
    R[1,1] =  c 
    
    return R 
    
def check_cone():
    
    
    
    a = 1.0
    r = 1.606587877768772
    N = 16
    
    h = 3.0
    total_posterior     = np.pi
    min_angle = -total_posterior/2;
    max_angle =  total_posterior/2;
    
    cone = cone_filter(a, r, h, N, min_angle, max_angle)
     
    X_flat = np.array([1.2, 1]) 
    
    print 'orig = ', X_flat
    
    X = cone.apply(X_flat) 
    
    print 'image = ', X
    
    X_flat_calc = cone.apply_inverse(X)
    
    print 'inverted = ', X_flat_calc
    
    
    '''
    N = 32; 
    
    cone_only = false; 
    if cone_only
    
        x = linspace(0,pi/2,N); 
    
        ring_quarter = [r*cos(x); r*sin(x); h*ones(size(x))]; 
        ring_plane = zeros(2,N); 
    
        for j=1:N
            ring_plane(:,j) = cone_filter_inv_right_cone(ring_quarter(:,j), filter_params); 
        end 
    
        fig = figure; 
        plot(ring_plane(1,:), ring_plane(2,:)); 
        axis equal 
        title('valve ring in plane, cone only')
    
        mesh = linspace(0,1,N); 
        rays_in_plane = zeros(2,N,N); 
    
        initial_pt = [a;0]; 
    
        for j=1:N
            for k=1:N
                rays_in_plane(:,j,k) = mesh(j)*initial_pt  + (1-mesh(j))*ring_plane(:,k); 
            end 
        end 
    
        fig = figure; 
        plot(ring_plane(1,:), ring_plane(2,:));
        hold on 
    
        for k=1:N
            plot(rays_in_plane(1,:,k), rays_in_plane(2,:,k)); 
        end 
        axis equal 
        title('rays in plane plus valve ring, cone only')
    
    
    
        rays_in_3d = zeros(3,N,N); 
        for j=1:N
            for k=1:N
                rays_in_3d(:,j,k) = cone_filter(rays_in_plane(1,j,k), rays_in_plane(2,j,k), filter_params); 
            end 
        end 
    
        fig = figure; 
        plot3(ring_quarter(1,:), ring_quarter(2,:), ring_quarter(3,:))
        hold on 
        for k=1:N
            plot3(rays_in_3d(1,:,k), rays_in_3d(2,:,k), rays_in_3d(3,:,k)); 
        end 
        axis equal
        title('cone part in space, cone only')
    
    end 
    
    
    tic 
    
    N = 32; 
    x = linspace(-3*pi/4,3*pi/4,N); 
    
    
    ring_half = [r*cos(x); r*sin(x); h*ones(size(x))]; 
    ring_plane = zeros(2,N); 
    
    for j=1:N
        ring_plane(:,j) = cone_filter_inv(ring_half(:,j), filter_params); 
    end 
    
    fig = figure; 
    plot(ring_plane(1,:), ring_plane(2,:)); 
    axis equal 
    title('valve ring in plane')
    
    mesh = linspace(0,1,N); 
    rays_in_plane_right = zeros(2,N,N); 
    rays_in_plane_left  = zeros(2,N,N); 
    
    initial_pt_right = [ a;0]; 
    initial_pt_left  = [-a;0]; 
    
    for j=1:N
        for k=1:N
            rays_in_plane_right(:,j,k) = mesh(j)*initial_pt_right  + (1-mesh(j))*ring_plane(:,k); 
            rays_in_plane_left(:,j,k)  = mesh(j)*initial_pt_left   + (1-mesh(j))*ring_plane(:,k); 
        end 
    end 
    
    fig = figure; 
    plot(ring_plane(1,:), ring_plane(2,:));
    hold on 
    
    for k=1:N
        plot(rays_in_plane_right(1,:,k), rays_in_plane_right(2,:,k)); 
    end 
    axis equal 
    title('rays in plane plus valve ring')
    
    
    
    rays_in_3d_right = zeros(3,N,N); 
    rays_in_3d_left  = zeros(3,N,N); 
    for j=1:N
        for k=1:N
            rays_in_3d_right(:,j,k) = cone_filter(rays_in_plane_right(1,j,k), rays_in_plane_right(2,j,k), filter_params); 
            rays_in_3d_left (:,j,k) = cone_filter(rays_in_plane_left (1,j,k), rays_in_plane_left (2,j,k), filter_params); 
        end 
    end 
    
    fig = figure; 
    plot3(ring_half(1,:), ring_half(2,:), ring_half(3,:))
    hold on 
    for k=1:N
        plot3(rays_in_3d_right(1,:,k), rays_in_3d_right(2,:,k), rays_in_3d_right(3,:,k)); 
        plot3(rays_in_3d_left (1,:,k), rays_in_3d_left (2,:,k), rays_in_3d_left (3,:,k)); 
    end 
    axis equal
    title('whole surface in space')
    
    
    
    % invert and make sure that we get something reasonable 
    rays_in_plane_right_computed = zeros(2,N,N); 
    rays_in_plane_left__computed = zeros(2,N,N); 
    for j=1:N
        for k=1:N
            rays_in_plane_right_computed(:,j,k) = cone_filter_inv(rays_in_3d_right(:,j,k), filter_params);  
            rays_in_plane_left__computed(:,j,k) = cone_filter_inv(rays_in_3d_left (:,j,k), filter_params);  
        end 
    end 
    
    fig = figure; 
    plot(ring_plane(1,:), ring_plane(2,:));
    hold on 
    
    for k=1:N
        plot(rays_in_plane_right_computed(1,:,k), rays_in_plane_right_computed(2,:,k)); 
        plot(rays_in_plane_left__computed(1,:,k), rays_in_plane_left__computed(2,:,k)); 
    end 
    axis equal 
    title('rays in plane plus valve ring - Computed by inversion')
    
    
    'total time for surf build'
    toc
    

    ''' 

 
