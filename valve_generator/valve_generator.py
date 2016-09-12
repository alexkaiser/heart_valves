'''
Generates initial mitral valve geometry in various configurations.


Alexander D. Kaiser
9/2016
'''


import numpy as np
import scipy
import reference_surface_generator
import math
from reference_surface_generator import cone_filter

class mitral_valve:

    def __init__(self, N):
        pass


class leaflet_and_chordeae:
    '''
    Data for leaflet and chordeae tendineae
    
    X Lea
    '''

    def __init__(self):
        self.X = None
        self.C = None
        self.R = None
        self.N = None
        self.left_free_edge  = None
        self.right_free_edge = None
        self.radial_and_circumfrential_fibers = False
        self.left_papillary  = None
        self.right_papillary = None



class leaflet:
    '''
    Leaflet class. Includes 
    '''
    
    
    def __init__(self, filter, radial_and_circumfrential_fibers = False):
        
        self.filter = filter 
        self.radial_and_circumfrential_fibers = False
        self.points = []
        self.N = filter.N


    def build_reference_surface_diag_fibers(self):
        '''
        Builds the reference surface associated with the current filter 
        '''
        
        assert not self.radial_and_circumfrential_fibers
        
        X      = np.zeros([3, self.N+1, self.N+1])
        X_flat = np.zeros([2, self.N+1, self.N+1])
        
        mesh_ring = np.linspace(self.filter.min_angle, self.filter.max_angle, self.N+1, endpoint=True) 
        
        ring      = np.zeros([3,self.N+1])
        ring[0,:] = r * np.cos(mesh_ring)
        ring[1,:] = r * np.sin(mesh_ring)
        ring[2,:] = self.filter.h * np.ones(N+1) 
        
        k = self.N
        for j in range(self.N+1): 
            X_flat[:,j,k] = self.filter.apply_inverse(ring[:,j])
            k -= 1
        
        
        
        
        
        
        
class chordeae_tree:

    def __init__(self):
        pass

class free_edge:
    ''' 
    Coordinates and indices of the free edge
    Used for connecting chordeae tendineae
    '''
    def __init__(self):
        pass

class counter:
    ''' 
    Counters on global and local indices
    '''
    def __init__(self, global_idx):
        self.global_idx = 0
        
    def intcrement(self):
        self.global_idx += 1 


class point_3:
    '''
    Three dimensional point class. 
    Includes 
    
    X             Coordinates 
    R             Refernce coordinate, if applicable 
    global_idx    Global index of point 
    local_idx     Local index, may be integer or tuple (in multi dimensional case)
    is_bc         Not included in updates to leaflet structure 
    is_leaflet    Is internal to leaflet 
    is_chordae    In chordae tree 
    '''
    
    def __init__(self):
        self.X          = None
        self.R          = None 
        self.global_idx = None
        self.local_idx  = None
        self.nbrs       = None 
        self.is_bc      = False 
        self.is_leaflet = False
        self.is_chordae = False




if __name__ == '__main__':

    print 'it is crowded...'


    pi = math.pi

    # filter constants
    a = 1.0
    r = 1.606587877768772
    N = 16
    
    h_posterior = 3.0
    total_posterior     = pi + pi/6 + pi/12;
    min_angle_posterior = -total_posterior/2;
    max_angle_posterior =  total_posterior/2;
    
    radial_and_circumfrential_fibers_posterior = False 
    
    filter_posterior = cone_filter(a, r, h_posterior, N, min_angle_posterior, max_angle_posterior)
    
    leaflet_posterior = leaflet(filter_posterior, radial_and_circumfrential_fibers_posterior)
    
    leaflet_posterior.build_reference_surface_diag_fibers()
    
    
    
    














