

'''
Takes a series of silo files containing Lagrangian mesh. 
Converts them into a single lines3d file.
Topology is read from spring file. 

Names are currently hardcoded, should fix this. 

'''

import os 



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
        
        # first two tokens are always the indices 
        indices = (int(line_split[0]), int(line_split[1]))
    
        spring_list.append(indices)

    return spring_list


def print_lines3d_header(lines3d_file, spring_list): 
    '''
    Prints the header for the list of springs 
    
    Particles not included yet   
    ''' 

    # always a 2 since springs have two ends  
    spring_str = '2 ' + str(len(spring_list)) + '\n'
    lines3d_file.write( spring_str )


def prepend_header(lines3d_file_name, spring_list, n_frames): 
    ''' 
    Adds the number of vertices and the number of frames to 
    the lines3d file. 
    
    Wasteful and out of place.
    
    Input: 
    lines3d_file    lines3d file without opening line 
    n_vertices      Number of vertices per frame 
    n_frames        Number of frames 
    ''' 

    temp_name = 'temp_file.txt' 
    temp = open(temp_name, 'w')
    lines_file = open(lines3d_file_name, 'r')

    # write the header 
    header = str(len(spring_list)) + ' ' + str(n_frames) + '\n'
    temp.write(header)
    
    # write the whole old file onto the new one 
    temp.write(lines_file.read())
    
    lines_file.close()
    temp.close()
    
    # clobber the old with the temp 
    os.rename(temp_name, lines3d_file_name)
    


def xyz_to_string_array(file, lines3d_file, spring_list, num_vertices_prev=None):
    '''
    Reads the vertices into a     
    
    Input: 
        file                File to read, must be in XYZ format 
        lines3d_file        Write to this file, header not included 
        num_vertices_prev   If not None, number reported in the file 
                            is required to equal this. 
    
    '''
    vertex_strings = []
    
    n_vertices = int(file.readline()) 
    
    # second line is always formatting junk
    file.readline()
    
    # check that two numbers of vertices are equal 
    if num_vertices_prev is not None: 
        assert n_vertices == num_vertices_prev 
            
    for line in file: 
        split = line.split()
        
        # always tokes 1,2,3, zero is a question mark 
        temp = split[1] + ' ' + split[2] + ' ' + split[3] + '\n'
        vertex_strings.append(temp)

    return vertex_strings, n_vertices
    
    
def write_vertices(lines3d_file, vertex_strings, spring_list): 
    '''
    Writes vertex strings according to what is listed 
    in the spring list. 
    
    Input: 
        lines3d_file       File Must be open. 
        vertex_strings     Vertices in string format
        spring_list        List of tuples of spring indices 
    '''
    
    for pair in spring_list: 
        lines3d_file.write(vertex_strings[pair[0]])
        lines3d_file.write(vertex_strings[pair[1]])




if __name__ == '__main__':

    print 'it is crowded...'

    spring_name = '../mitral_tree.spring'
    spring_list = read_springs(spring_name)

    lines3d_file_name = 'mitral_tree.3D'
    lines3d_file      = open(lines3d_file_name, 'w')
    
    # print the header 
    print_lines3d_header(lines3d_file, spring_list)

    n_vertices = None
    n_frames = 0 

    # for all the sorted files... 
    for f_name in sorted(os.listdir(os.getcwd())): 
        if f_name.endswith('.xyz'):
                        
            # have a valid file 
            n_frames += 1 
    
            temp_file = open(f_name, 'r')
            
            vertex_strings, n_vertices = xyz_to_string_array(temp_file, lines3d_file, spring_list, n_vertices)
            write_vertices(lines3d_file, vertex_strings, spring_list) 
            
            temp_file.close()
            
            
    
    # close before messing with things  
    lines3d_file.close()
    
    # fix the header with the number of fibers (here all length two springs) and frames 
    prepend_header(lines3d_file_name, spring_list, n_frames)

    print 'but we out here'






