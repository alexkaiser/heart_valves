

'''
Takes a series of silo files containing Lagrangian mesh. 
Converts them into a single lines3d file.
Topology is read from spring file. 

'''

# Copyright (c) 2019, Alexander D. Kaiser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os 
import sys 
import numpy as np
from multiprocessing import Process

def read_springs(spring_name, target_list=None): 
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
        indices = (int(line_split[0]), int(line_split[1]), float(line_split[2]), float(line_split[3]) , int(line_split[4]))

        spring_list.append(indices)

    if target_list is not None: 
        spring_list_cropped = []

        for spring_data in spring_list:
            idx          = spring_data[0]
            idx_nbr      = spring_data[1]
            function_idx = spring_data[4]

            # output vertices that are targets 
            # and vertices that have function_idx == 1
            # which is 
            if function_idx == 1:
                if (idx in target_list) or (idx_nbr in target_list):

                    if (idx in target_list) and (idx_nbr in target_list):
                        raise ValueError('both idx and nbr should not be targets')

                    spring_list_cropped.append(spring_data)

        sprint_list = spring_list_cropped

    return spring_list


def read_targets(target_name):

    target = open(target_name, 'r')

    # first line is always the number of springs 
    n_targets = int(target.readline()) 

    target_list = []

    for i in range(n_targets): 
        line_split = (target.readline()).split()
        idx = int(line_split[0])
        target_list.append(idx)

    return target_list


def xyz_to_string_array(file, spring_list, num_vertices_prev=None):
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
    
    
def write_vertices(lines3d_file, vertex_strings, spring_list, n_particles=None): 
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
        
    if n_particles is not None: 
        for i in range(len(vertex_strings) - n_particles, len(vertex_strings)): 
            lines3d_file.write(vertex_strings[i])


def compute_forces_and_tangents(output_file, spring_list, target_list, vertex_strings):

    for spring_data in spring_list:
        idx          = spring_data[0]
        idx_nbr      = spring_data[1]
        kappa        = spring_data[2]
        rest_len     = spring_data[3]
        function_idx = spring_data[4]

        found_valid = False 

        # output vertices that are targets 
        # and vertices that have function_idx == 1
        # which is 
        if function_idx == 1:
            if (idx in target_list) or (idx_nbr in target_list):

                if (idx in target_list) and (idx_nbr in target_list):
                    raise ValueError('both idx and nbr should not be targets')

                if idx in target_list:
                    target_idx = idx
                    pair_idx = idx_nbr

                elif idx_nbr in target_list:
                    target_idx = idx_nbr
                    pair_idx = idx

                coords = (vertex_strings[target_idx]).split()
                position_target = np.array([float(c) for c in coords])

                coords_pair = (vertex_strings[pair_idx]).split()
                position_pair = np.array([float(c) for c in coords_pair])

                # print 'position_target = ', position_target
                # print 'position_pair = ', position_pair

                tangent = (position_pair - position_target) / np.linalg.norm(position_pair - position_target) 

                force_magnitude = spring_function_collagen(position_target, position_pair, kappa, rest_len)

                force = force_magnitude * tangent 

                to_write  = str(target_idx) + ' ' + str(pair_idx) + ' '
                for c in coords:
                    to_write += c + ' '
                for f in force:
                    to_write += str(f) + ' '
                to_write += '\n'

                output_file.write(to_write)


def spring_function_collagen(position, position_nbr, kappa, rest_len):
    '''
    Compute force for collagen springs
    
    Zero force under compression
    
    Exponential growth until full recruitment,
    experimentally determined strain at which microscopic fibers align
     
    Linear after full recruitment
    '''

    MPa_TO_CGS           = 1.0e7    
    a                    = 4643.4       # Coeff of exponential term
    b                    = 49.9643      # Exponential rate
    full_recruitment     = 0.145             # Linear at strains larger than this
    eta_collagen         = 32.5 * MPa_TO_CGS # Linear slope, in barye = dynes/cm^2 = g cm/(s cm^2)
    collagen_x_intercept = 0.125             # Linear collagen part intercepts x axis at this strain
    collagen_y_intercept = -collagen_x_intercept * eta_collagen # Linear collagen part intercepts y axis at this stress

    R = np.linalg.norm(position - position_nbr)
    
    # Strain, dimension 1
    E = R/rest_len - 1.0
    
    # Compute the force
    if (E > full_recruitment):
        return kappa * (eta_collagen*E + collagen_y_intercept)
    elif (E > 0.0):
        return kappa * a * (np.exp(b*E) - 1.0)
    else:
        return 0.0

    return 0.0




if __name__ == '__main__':

    print 'it is crowded...'

    if len(sys.argv) <= 1:
        print 'defaulting to default file name'
        base_name = "mitral_tree"
    else: 
        base_name = str(sys.argv[1])
        
    if len(sys.argv) >= 4:
        total_procs = int(sys.argv[2])
        proc_num    = int(sys.argv[3])
    else: 
        total_procs = 1
        proc_num = 0 
    
    spring_name = base_name + '.spring'
    spring_list = read_springs(spring_name)

    target_name = base_name + '.target'
    target_list = read_targets(target_name)

    n_vertices = None
    valid_loop_count = 0

    # for all the sorted files... 
    for f_name in sorted(os.listdir(os.getcwd())): 
        if f_name.startswith(base_name) and f_name.endswith('.xyz'):

            if (valid_loop_count % total_procs) == proc_num:

                print 'printing frame', str(valid_loop_count), '\n'
                        
                with open(f_name, 'r') as f:
                    vertex_strings, n_vertices = xyz_to_string_array(f, spring_list, n_vertices)

                frame_num_as_str = '{0:05d}'.format(valid_loop_count)
                file_name = 'force_data_' + frame_num_as_str + '.txt'

                with open(file_name, 'w') as f: 
                    compute_forces_and_tangents(f, spring_list, target_list, vertex_strings)

            valid_loop_count += 1 
    
    print 'but we out here'

