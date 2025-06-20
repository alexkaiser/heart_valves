
import os 
import re 
import numpy as np
import sys 
import pyvista 

def process_vtk_file_pyvista(file_name, base_name_out, scaling=1.0):


    vertex_file = open(base_name_out + '.vertex', 'w')

    mesh = pyvista.read(file_name)

    vertex_file.write(str(mesh.n_points) + '\n')

    vertices = scaling * mesh.points 

    print( 'vertices = ', vertices[0][:])
    print( 'vertices = ', vertices[1][:])

    spring_idx = []
    n_springs = 0
    first_placed = None 
    n_placed = 0

    for i in range(mesh.n_cells):
        if mesh.cells[3*i] != 2: 
                raise ValueError('Must have edges (not triangles or other data structures)')
        spring_idx.append( (mesh.cells[3*i + 1], mesh.cells[3*i + 2]) )


    current_idx = spring_idx[0][0]
    first_idx   = current_idx
    next_idx    = spring_idx[0][1]
    complete    = False; 
    while True: 

        for d in range(3):
            vertex_file.write(str(vertices[current_idx][d]) + ' ')
        vertex_file.write('\n')

        found_next = False 
        for pair in spring_idx:
            if pair[0] == next_idx:
                found_next  = True 
                current_idx = next_idx
                next_idx    = pair[1]
                break 

        if not found_next:
            raise ValueError('Did not find a vertex with the correct index')

        if current_idx == first_idx:
            complete = True
            break 

    vertex_file.close()


def process_vtk_file(file_name, base_name_out, scaling=1.0):
    f = file_name

    strings = f.split('.')
    print('found vtk file ', f )

    vtk_file  = open(f, 'r')
    vtk_file_as_string = vtk_file.read()
    vtk_file.close()

    res = re.search('POINTS(.*)float', vtk_file_as_string)
    n_points = int(res.group(1))
    print( "n_points = ", n_points)

    vertex_file = open(base_name_out + '.vertex', 'w')

    # header ends with float\n
    # header into vertices_split_1[0]
    # rest of file into vertices_split_1[1]
    vertices_split_1 = vtk_file_as_string.split('float')
    # next thing after vertices is 'METADATA'
    if 'METADATA' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('METADATA')
    elif 'CELLS' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('CELLS')
    else:
        raise ValueError("Vertices must end with METADATA or CELLS and otherwise parsing not implemented")

    all_vertices_string = vertices_split_2[0]

    # may or may not be leading white space, remove if so 
    all_vertices_string = all_vertices_string.lstrip()

    vertex_file.write(str(n_points) + '\n')

    vertices = []
    i = 0
    for coord in all_vertices_string.split():
        # coords to floats for spring files 
        vertices.append(scaling * float(coord))

        # # write exact string to file 
        # vertex_file.write(coord + " ")
        # i += 1 
        # if (i%3) is 0:
        #     vertex_file.write("\n")

    vertices = np.reshape(vertices, (n_points,3))

    print( 'vertices = ', vertices[0][:])
    print( 'vertices = ', vertices[1][:])

    '''    
    for pt in range(n_points):
        for d in range(3):
            vertex_file.write(str(vertices[pt][d]) + ' ')
        vertex_file.write('\n')

    vertex_file.close()
    '''

    # just reopen becuase why not 
    vtk_file  = open(f, 'r')
    writing = False
    spring_idx = []
    n_springs = 0
    first_placed = None 
    n_placed = 0

    for line in vtk_file: 
        if writing:

            # print 'line.split = ', line.split()

            # line has to start with 4 or end writing 
            # assuming tetrahedral mesh 
            # print 'line = ', line 
            line_split = line.split()
            
            # break on empty line, or whitespace only line 
            if len(line_split) == 0:
                break

            # also break if first token is not 4
            first_token = line_split[0]
            # print 'first_token = ', first_token

            try:
                n_edges = int(first_token)
            except ValueError:
                print( "Not an int, break")
                break 

            if first_token == "CELL_TYPES": 
                break

            if n_edges != 2: 
                raise ValueError('Must have edges (not triangles or other data structures)')

            spring_idx.append( (int(line_split[1]), int(line_split[2])))
            n_springs += 1
            
        if line.startswith('CELLS') or line.startswith('LINES'):
            writing = True 
            lines = 0


    current_idx = spring_idx[0][0]
    first_idx   = current_idx
    next_idx    = spring_idx[0][1]
    complete    = False; 
    while True: 

        for d in range(3):
            vertex_file.write(str(vertices[current_idx][d]) + ' ')
        vertex_file.write('\n')

        found_next = False 
        for pair in spring_idx:
            if pair[0] == next_idx:
                found_next  = True 
                current_idx = next_idx
                next_idx    = pair[1]
                break 

        if not found_next:
            raise ValueError('Did not find a vertex with the correct index')

        if current_idx == first_idx:
            complete = True
            break 

    vertex_file.close()
    vtk_file.close()


if __name__== "__main__":

    if len(sys.argv) >= 2:
        file_name = str(sys.argv[1])

        if len(sys.argv) >= 3:
            scaling = float(sys.argv[2])
        else:
            scaling = 1.0 

        assert file_name.endswith('.vtk')
        base_name_out = file_name.rsplit('.',1)[0]
        process_vtk_file(file_name, base_name_out, scaling)


    # further scripting if desired

    inlet_outlet_files = False 
    if inlet_outlet_files:
        scaling = .1            
        # inlet outlet files are vertex files only 
        spring_strength_rel = 0.0
        target_strength  = 0.0
        damping_strength = 0.0
        
        file_name = "aorta_bdry.vtk"
        base_name_out = "aorta_bdry"
        process_vtk_file(file_name, base_name_out, scaling)            

        file_name = "atrium_bdry.vtk"
        base_name_out = "atrium_bdry"
        process_vtk_file(file_name, base_name_out, scaling)            

        # file_name = "aortic_annulus.vtk"
        # base_name_out = "aortic_annulus"
        # process_vtk_file(file_name, base_name_out, scaling)            

    mri4d_compare_aortic_files = False  
    if mri4d_compare_aortic_files:
        # inlet outlet files are vertex files only 
        scaling = 1
        
        # file_name = "right_ventricle.vtk"
        # base_name_out = "right_ventricle_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)            

        # file_name = "right_pa.vtk"
        # base_name_out = "right_pa_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)

        # file_name = "left_pa.vtk"
        # base_name_out = "left_pa_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)

        # file_name = "right_ventricle_pt025cm_inwards.vtk"
        # base_name_out = "right_ventricle_pt025cm_inwards_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)            

        # file_name = "right_ventricle_pt025cm_out.vtk"
        # base_name_out = "right_ventricle_pt025cm_out_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)            

        # file_name = "right_ventricle_pt05cm_out.vtk"
        # base_name_out = "right_ventricle_pt05cm_out_bdry"
        # process_vtk_file(file_name, base_name_out, scaling)            

        file_name = "right_ventricle_pt1cm_out.vtk"
        base_name_out = "right_ventricle_pt1cm_out_bdry"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)  

    new_parser_test = False 
    if new_parser_test:
        scaling = 1 

        file_name = "left_pa.vtk"
        base_name_out = "left_pa_bdry_orig"
        process_vtk_file(file_name, base_name_out, scaling)

        file_name = "left_pa.vtk"
        base_name_out = "left_pa_bdry_pyvista"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

    aortic_384_bicuspidization_files = False 
    if aortic_384_bicuspidization_files:
        scaling = 0.1

        file_name = 'aorta_bdry_384_layer_3_constriction.vtu'
        base_name_out = "aorta_bdry_384"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

        file_name = 'lvot_bdry_384_layer_3_constriction.vtu'
        base_name_out = "lvot_bdry_384"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

    aortic_192_bicuspidization_files = False 
    if aortic_192_bicuspidization_files:
        scaling = 0.1

        file_name = "aorta_bdry_192_layer_3_double_constriction.vtu"
        base_name_out = "aorta_bdry_192"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

        file_name = "lvot_bdry_192_layer_3_constriction.vtu"
        base_name_out = "lvot_bdry_192"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

    seg_test = True 
    if seg_test:
        scaling = 0.1

        file_name = 'aorta_bdry_outlet.vtu'
        base_name_out = "aorta_bdry_384"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

        file_name = 'lvot_bdry_constriction.vtu'
        base_name_out = "lvot_bdry_384"
        process_vtk_file_pyvista(file_name, base_name_out, scaling)

