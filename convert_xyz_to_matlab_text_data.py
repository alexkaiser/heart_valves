

'''
Takes a series of silo files containing Lagrangian mesh. 
Converts them into a single lines3d file.
Topology is read from spring file. 

'''

import os 
import sys 
import math 

DEBUG = True


def read_springs(spring_name): 
    ''' 
    Reads a file in IBAMR .spring format. 
    Returns a list of length 2 tuples, 
    each of which are the indices of a spring in the simulation.
    
    Input: 
        spring_name    File name, must be in current directory 
       
    Output: 
        array of length two tuples of ints 
    '''

    spring = open(spring_name, 'r')

    # first line is always the number of springs 
    n_springs = int(spring.readline()) 

    spring_list = []

    
    for i in range(n_springs): 
        line_split = (spring.readline()).split()
        
        include = True 
        
        # first two tokens are always the indices 
        indices = (int(line_split[0]), int(line_split[1]))

        # token after comment may contain information on whether to include this spring 
        if '#' in line_split:
            idx = line_split.index('#') + 1
            if idx < len(line_split):
                include = int(line_split[idx])
            
        # include = True
        if include: 
            spring_list.append(indices)

    return spring_list
    


def xyz_to_vertices_list(file, spring_list, num_vertices_prev=None):
    '''
    Reads the vertices into an array of tupeles      
    
    Input: 
        file                File to read, must be in XYZ format 
        lines3d_file        Write to this file, header not included 
        num_vertices_prev   If not None, number reported in the file 
                            is required to equal this. 
    
    '''
    vertex_coords = []
    
    n_vertices = int(file.readline()) 
    
    # second line is always formatting junk
    file.readline()
    
    # check that two numbers of vertices are equal 
    if num_vertices_prev is not None: 
        assert n_vertices == num_vertices_prev 
            
    for line in file: 
        split = line.split()
        
        # always tokens 1,2,3, token zero is a question mark 
        temp = (split[1], split[2], split[3])
        
        vertex_coords.append(temp)

    # print 'vertex_coords = ', vertex_coords
    # print 'len(vertex_coords) = ', len(vertex_coords)

    return vertex_coords, n_vertices


def xyz_to_particle_vertices_only(file, n_vertices, n_particles):
    '''
    Reads the vertices into an array of tupeles      
    
    Input: 
        file                File to read, must be in XYZ format 
        lines3d_file        Write to this file, header not included 
        num_vertices_prev   If not None, number reported in the file 
                            is required to equal this. 
    
    '''
    vertex_coords = []
    
    # n_vertices = int(file.readline())     
    # second line is always formatting junk
    # file.readline()
            
    min_particle_idx = n_vertices - n_particles
    min_particle_line_number = min_particle_idx + 2 # two header lines in file that are ignored

    n_particles_placed = 0

    for i, line in enumerate(file):

        if i >= min_particle_line_number:
            split = line.split()
            
            # always tokens 1,2,3, token zero is a question mark 
            temp = (split[1], split[2], split[3])
            
            vertex_coords.append(temp)

            n_particles_placed += 1 

    assert n_particles == n_particles_placed

    return vertex_coords



    
def write_spring_vertices(matlab_data_file, vertex_coords, spring_list): 
    '''
    Writes vertex strings according to what is listed 
    in the spring list. 
    
    Input: 
        lines3d_file       File Must be open. 
        vertex_strings     Vertices in string format
        spring_list        List of tuples of spring indices 
    '''
    
    names = ['x_coords_springs', 'y_coords_springs', 'z_coords_springs']
    dimension = 3 

    for dim in range(dimension):
        matlab_data_file.write(names[dim] + ' = [')

        for pair in spring_list: 

            #print 'pair[0] = ', pair[0]
            #print 'pair[1] = ', pair[1]
            #print 'len(vertex_coords) = ', len(vertex_coords)

            string_lower = vertex_coords[pair[0]][dim]
            string_upper = vertex_coords[pair[1]][dim]
            matlab_data_file.write(string_lower + ' ' + string_upper + ';\n')

        matlab_data_file.write(']\';\n\n')

    '''        
    if n_particles is not None: 
        for i in range(len(vertex_strings) - n_particles, len(vertex_strings)): 
            lines3d_file.write(vertex_strings[i])
    '''



def get_full_xyz_name(xyz_base_name, frame_number, data_dir='.', xyz_fill_len=10):
    return data_dir + '/' + xyz_base_name + str(frame_number).zfill(xyz_fill_len) + '.xyz'


def write_comet_tail_data(matlab_data_file, frame_number, n_vertices, n_particles, xyz_base_name, comet_tail_len, dt, min_frame, n_steps): 
    '''
    Writes vertex strings according to what is listed 
    in the spring list. 
    
    Input: 
        lines3d_file       File Must be open. 
        vertex_strings     Vertices in string format
        spring_list        List of tuples of spring indices 
    '''
    
    #matlab_data_file.write('comet_tail_len = ' + str(comet_tail_len) + ';\n')
    #matlab_data_file.write('x_comet_coords = zeros(n_particles, comet_tail_len)'); 
    #matlab_data_file.write('y_comet_coords = zeros(n_particles, comet_tail_len)'); 
    #matlab_data_file.write('z_comet_coords = zeros(n_particles, comet_tail_len)'); 
    
    '''
    if n_particles is not None: 
        for i in range(len(vertex_strings) - n_particles, len(vertex_strings)): 
            lines3d_file.write(vertex_strings[i])
    '''

    particle_vertices = []
    dimension = 3

    nans = []
    for i in range(n_particles):
        nans.append( ('NaN', 'NaN', 'NaN'))

    for frame in range(frame_number, frame_number - comet_tail_len, -1):
        if min_frame <= frame < n_steps: 
            xyz_file_name = get_full_xyz_name(xyz_base_name, frame)
            with open(xyz_file_name) as xyz_file: 
                particle_vertices.append(xyz_to_particle_vertices_only(xyz_file, n_vertices, n_particles))
        else:
            # particles are out of range of the simulation 
            particle_vertices.append(nans)

    # for centered velocity at PRESENT time step     
    if (frame_number+1) < n_steps:
        future_file_name = get_full_xyz_name(xyz_base_name, frame_number+1)
        with open(future_file_name) as xyz_file:
            future_coords = xyz_to_particle_vertices_only(xyz_file, n_vertices, n_particles)
    else: 
        assert False 
        future_coords = nans

    # for centered velocity at the end of the comet tail 
    if (frame_number - comet_tail_len) > 0:
        past_file_name = get_full_xyz_name(xyz_base_name, frame_number - comet_tail_len)
        with open(past_file_name) as xyz_file:
            past_coords = xyz_to_particle_vertices_only(xyz_file, n_vertices, n_particles)
    else: 
        # assert False 
        past_coords = nans


    velocity_norms = [['' for particle_idx in range(n_particles)] for tail_idx in range(comet_tail_len)]

    # if any link is greater than this in norm, the remainder of the tail is replaced with NaNs 
    norm_threshold = 5
    for particle_idx in range(n_particles):
        tail_valid = True 
        # compare current to prevoius starting at second (index 1)
        # set current coords to NaN if not 
        for tail_idx in range(0,comet_tail_len):

            norm = 0.0
            velocity = 0.0
            for dim in range(dimension):

                if tail_idx == (comet_tail_len-1):
                    next_pos = float(past_coords[particle_idx][dim])
                else:
                    next_pos = float(particle_vertices[tail_idx+1][particle_idx][dim])

                val          = float(particle_vertices[tail_idx  ][particle_idx][dim])

                if tail_idx == 0:
                    prev     = float(future_coords[particle_idx][dim])
                else:
                    prev     = float(particle_vertices[tail_idx-1][particle_idx][dim])

                # for cropping periodically use absolute distance 
                norm  += (val - next_pos)**2

                # velocity norm gets a centered difference 
                velocity += (next_pos - prev)**2
            
            norm = math.sqrt(norm)
            if norm >= norm_threshold:
                tail_valid = False

            velocity = math.sqrt(velocity) / (2.0 * dt)
            velocity_norms[tail_idx][particle_idx] = str(velocity)

            if not tail_valid:
                particle_vertices[tail_idx][particle_idx] = ('NaN', 'NaN', 'NaN')
                velocity_norms[tail_idx][particle_idx] = 'NaN'


    # particle_vertices is a 2d list of tuples 
    # first idx is number of time steps from current point 
    # next idx is particle idx 
    # final is tuple index for coordinate 
    names = ['x_comet_coords', 'y_comet_coords', 'z_comet_coords', 'particle_velocity']

    # each dimension separately 
    for dim,name in enumerate(names):
        matlab_data_file.write(name + ' = [')
        for particle_idx in range(n_particles):
            for tail_idx in range(comet_tail_len):
                if name is 'particle_velocity':
                    matlab_data_file.write(velocity_norms[tail_idx][particle_idx] + ' ')                
                else:
                    matlab_data_file.write(particle_vertices[tail_idx][particle_idx][dim] + ' ')
            matlab_data_file.write(';\n')
        matlab_data_file.write(']\';\n\n')


if __name__ == '__main__':

    print 'it is crowded...'

    print 'Default arguments:'
    print 'base_name = "mitral_tree_512" '
    print 'frame_number = 100'
    print 'comet_tail_len = 20'
    print 'dt = 1.5e-6'

    base_name = 'mitral_tree_512'
    xyz_base_name = 'mitral_tree_512_lines3d_'
    xyz_fill_len = 10
    # frame_number = 1005
    frames = [33] # range(1441) 
    # frames = [1005, 1267, 1382]
    comet_tail_len = 20
    dt_sim = 1.5e-6
    output_frequency = 1111 
    dt = dt_sim * output_frequency
    data_dir = '.'
    min_step = 0
    n_steps  = 1442


    spring_name = base_name + '.spring'
    spring_list = read_springs(spring_name)
        
    # find out how many particles there are
    # particles are always placed last 
    try: 
        particles_file_name = base_name + '.particles'
        particles_file = open(particles_file_name, 'r')
        n_particles = int(particles_file.readline())
        particles_file.close()
    except: 
        n_particles = None

    for frame_number in frames:

        matlab_data_file_name  = base_name + '_' + str(frame_number).zfill(6) + '.m'
        matlab_data_file       = open(matlab_data_file_name, 'w')

        n_vertices = None
        n_frames = 0 
        valid_loop_count = 0

        '''
        # for all the sorted files... 
        for f_name in sorted(os.listdir(os.getcwd())): 
            if f_name.startswith(base_name) and f_name.endswith('.xyz'):
                if (valid_loop_count % frame_stride) == 0:                
                    print 'printing frame', str(valid_loop_count), '\n'
                            
                    # have a valid file 
                    n_frames += 1 
        
                    temp_file = open(f_name, 'r')
                
                    vertex_strings, n_vertices = xyz_to_string_array(temp_file, lines3d_file, spring_list, n_vertices)
                    write_vertices(lines3d_file, vertex_strings, spring_list, n_particles)             
                    temp_file.close()
                valid_loop_count += 1 
        '''

        xyz_file_name = get_full_xyz_name(xyz_base_name, frame_number, data_dir, xyz_fill_len)
        with open(xyz_file_name) as xyz_file:
            vertex_coords, n_vertices = xyz_to_vertices_list(xyz_file, spring_list)
            # print 'vertex_coords = ', vertex_coords
            print 'len(vertex_coords) = ', len(vertex_coords)
            print 'n_vertices = ', n_vertices
            write_spring_vertices(matlab_data_file, vertex_coords, spring_list) 


        write_comet_tail_data(matlab_data_file, frame_number, n_vertices, n_particles, xyz_base_name, comet_tail_len, dt, min_step, n_steps)
     
        matlab_data_file.close()

    print 'but we out here'






