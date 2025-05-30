
import os 
import re 
import numpy as np
import pyvista


def spring_string(indices, vertices, spring_strength_rel, use_k_abs=False):

    tol = 1e-10
    v0 = np.array(vertices[indices[0]][:])
    v1 = np.array(vertices[indices[1]][:])
    rest_len = np.linalg.norm(v0 - v1)
    if (spring_strength_rel > 0.0) and (rest_len > tol):
        k_abs = spring_strength_rel / rest_len
    else:
        k_abs = 0.0

    # print("v0 = ", v0, "v1 = ", v1, "k_abs = ", k_abs, "rest_len = ", rest_len)

    return str(indices[0]) + ' ' + str(indices[1]) + ' ' + str(k_abs) + ' ' + str(rest_len) + '\n'
    #return f'{indices[0]:d} {indices[1]:d} {k_abs:.10f} {rest_len:.10f} \n'

def get_hex_pair_list():

    # horizontal 
    pairs = []
    for j in range(4):
        pairs.append( (j%4, (j+1)%4) )
    for j in range(4):
        pairs.append( ((j%4) + 4, ((j+1)%4) + 4) ) 

    # vertical 
    for j in range(4):
        pairs.append((j,j+4))

    # four diagonals 
    pairs.append((0,6))    
    pairs.append((1,7))    
    pairs.append((2,4))    
    pairs.append((3,5))

    return pairs


# def get_cells(mesh):
#     """Returns a list of the cells from this mesh.
#     This properly unpacks the VTK cells array.
#     There are many ways to do this, but this is
#     safe when dealing with mixed cell types."""

#     # source: 
#     # https://github.com/pyvista/pyvista-support/issues/66

#     mesh_unstructured = pyvista.UnstructuredGrid(mesh)

#     offset = 0
#     cells = []
#     for i in range(mesh_unstructured.n_cells):
#         loc = i + offset
#         nc = mesh_unstructured.cells[loc]
#         offset += nc
#         cell = mesh_unstructured.cells[loc+1:loc+nc+1]
#         cells.append(cell)
#     return cells


def process_vtk_file_pyvista(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength=0.0, scaling=1.0, zero_springs=False):


    mesh = pyvista.read(file_name)

    mesh_unst = pyvista.UnstructuredGrid(mesh)
    mesh_unst.save(base_name_out + ".vtu")

    f = file_name

    vertex_file = open(base_name_out + '.vertex', 'w')

    vertex_file.write(str(mesh.n_points) + '\n')

    if (spring_strength_rel > 0.0) or zero_springs:
        spring_file = open(base_name_out + '.spring', 'w')
    
    if target_strength > 0.0:
        target_file = open(base_name_out + '.target', 'w')

    vertices = scaling * mesh.points 

    print( 'vertices = ', vertices[0][:])
    print( 'vertices = ', vertices[1][:])

    if target_strength > 0.0:
        target_file.write(str(mesh.n_points) + '\n')
    
    for pt in range(mesh.n_points):
        for d in range(3):
            vertex_file.write(str(vertices[pt][d]) + ' ')
        vertex_file.write('\n')

        if target_strength > 0.0:
            target_str = str(pt) + ' ' + str(target_strength)
            if damping_strength > 0.0:
                target_str += ' ' + str(damping_strength)
            target_file.write(target_str + '\n')

    vertex_file.close()

    if target_strength > 0.0:
        target_file.close()

    if (spring_strength_rel > 0.0) or zero_springs: 

        # just reopen becuase why not 
        writing = False
        spring_idx = set()
        n_springs = 0

        hex_pairs = get_hex_pair_list()

        # calling cells with a polydata or other type seems to fail 
        # but casting to unstructured first seems to fix 
        cells = pyvista.UnstructuredGrid(mesh).cells

        # print("cells = ", cells)

        offset = 0
        for i in range(mesh.n_cells):

            loc = i + offset
            nc = cells[loc]
            offset += nc
            this_cell_incides = cells[loc+1:loc+nc+1]

            # print("this_cell_incides = ", this_cell_incides)

            # hex mesh, do not do all to all 
            if len(this_cell_incides) == 8: 
                for pair in hex_pairs:
                    indices = (this_cell_incides[pair[0]], this_cell_incides[pair[1]])
                    if indices not in spring_idx:
                        spring_idx.add(indices)
                        n_springs += 1

            else:
                # default is all to all 
                # all pairs of indices (but only once)
                for i in this_cell_incides:
                    for j in this_cell_incides:
                        if i<j: 
                            if (i,j) not in spring_idx:
                                spring_idx.add((i,j))
                                n_springs += 1

        spring_file.write(str(n_springs) + '\n')
        sprint_strings = []
        for pair in spring_idx:
            if (spring_strength_rel > 0.0) or zero_springs:
                spring_file.write(spring_string(pair, vertices, spring_strength_rel))

        spring_file.close()




def process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength=0.0, scaling=1.0, zero_springs=False):
    f = file_name

    strings = f.split('.')
    print ('found vtk file ', f )

    vtk_file  = open(f, 'r')
    vtk_file_as_string = vtk_file.read()
    vtk_file.close()

    res = re.search('POINTS(.*)(float|double)', vtk_file_as_string)
    n_points = int(res.group(1))
    print( "n_points = ", n_points)

    vertex_file = open(base_name_out + '.vertex', 'w')
    
    if (spring_strength_rel > 0.0) or zero_springs:
        spring_file = open(base_name_out + '.spring', 'w')
    
    if target_strength > 0.0:
        target_file = open(base_name_out + '.target', 'w')

    # header ends with float\n
    # header into vertices_split_1[0]
    # rest of file into vertices_split_1[1]
    vertices_split_1 = re.split('float|double', vtk_file_as_string)
    # vtk_file_as_string.split('float')
    # next thing after vertices is 'METADATA'
    if 'METADATA' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('METADATA')
    elif 'CELLS' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('CELLS')
    elif 'VERTICES' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('VERTICES')
    elif 'POLYGONS' in vertices_split_1[1]:
        vertices_split_2 = (vertices_split_1[1]).split('POLYGONS')
    else:
        raise ValueError("Vertices must end with METADATA, CELLS, VERTICES or POLYGONS and otherwise parsing not implemented")

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

    if target_strength > 0.0:
        target_file.write(str(n_points) + '\n')
    
    for pt in range(n_points):
        for d in range(3):
            vertex_file.write(str(vertices[pt][d]) + ' ')
        vertex_file.write('\n')

        if target_strength > 0.0:
            target_str = str(pt) + ' ' + str(target_strength)
            if damping_strength > 0.0:
                target_str += ' ' + str(damping_strength)
            target_file.write(target_str + '\n')

    vertex_file.close()

    if target_strength > 0.0:
        target_file.close()

    if (spring_strength_rel > 0.0) or zero_springs: 

        # just reopen becuase why not 
        vtk_file  = open(f, 'r')
        writing = False
        spring_idx = set()
        n_springs = 0

        hex_pairs = get_hex_pair_list()

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
                    print ("Not an int, break, line = ", line )
                    break 

                if first_token == "CELL_TYPES": 
                    break

                this_cell_incides = []
                for j in range(1,n_edges+1):
                    this_cell_incides.append(int(line_split[j])) 

                # hex mesh, do not do all to all 
                if n_edges == 8: 
                    for pair in hex_pairs:
                        indices = (this_cell_incides[pair[0]], this_cell_incides[pair[1]])
                        if indices not in spring_idx:
                            spring_idx.add(indices)
                            n_springs += 1

                else:
                    # default is all to all 
                    # all pairs of indices (but only once)
                    for i in this_cell_incides:
                        for j in this_cell_incides:
                            if i<j: 
                                if (i,j) not in spring_idx:
                                    spring_idx.add((i,j))
                                    n_springs += 1
                
            if line.startswith('CELLS') or line.startswith('POLYGONS'):
                writing = True 
                lines = 0

        spring_file.write(str(n_springs) + '\n')
        sprint_strings = []
        for pair in spring_idx:
            if (spring_strength_rel > 0.0) or zero_springs:
                spring_file.write(spring_string(pair, vertices, spring_strength_rel))

        spring_file.close()

    vtk_file.close()


if __name__== "__main__":

    spring_strength_rel = 0.0

    # papillary target strength at 256 
    # this shoudn't change 
    target_strength  = 4.586778715596331e+05
    damping_strength = 0.0

    file_name = "ventricle.vtk"

    base_name = "LV_pt05cm_normal_1" #
    # base_name = "LV_pt05cm_hocm_d_"

    scaling = .1

    lv_files = False

    if lv_files:
        # initial file including target points etc 
        zero_springs = True
        process_vtk_file(file_name, base_name, spring_strength_rel, target_strength, damping_strength, scaling, zero_springs)

    registered_files = False
    if registered_files:

        total_frames = 20

        registered_base_names = ['ventricle_points_registered_masked_', 
                                 'mitral_no_partition_256_points_registered_masked_',
                                 'aortic_no_partition_192_points_registered_masked_']

        vertex_basename_out = ['LV_pt05cm_normal_1', 
                               'mitral_no_partition_256', 
                               'aortic_no_partition_192']

        for reg_base_name, output_base_name in zip(registered_base_names, vertex_basename_out):

            # registered files are vertex files only 
            spring_strength_rel = 0.0
            target_strength  = 0.0
            damping_strength = 0.0

            for number in range(total_frames):

                file_name = reg_base_name + str(number) + ".vtk"
                base_name_out = output_base_name + str(number).zfill(2) # no extension here 

                print ("file_name = ", file_name)
                print ("base_name_out = ", base_name_out)

                process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)            

    inlet_outlet_files = False 
    if inlet_outlet_files:
        # inlet outlet files are vertex files only 
        spring_strength_rel = 0.0
        target_strength  = 0.0
        damping_strength = 0.0
        
        file_name = "aorta_bdry.vtk"
        base_name_out = "aorta_bdry"
        process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)            

        file_name = "atrium_bdry.vtk"
        base_name_out = "atrium_bdry"
        process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)            

    mri4d_compare_aortic_files = False 
    if mri4d_compare_aortic_files:

        # these are already in CGS units 
        scaling = 1.0

        # target strength aortic_192
        target_strength = 388196.97181458189152

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.0147
        kappa_abs = 0.01 * 4.2649e+06 / 3.0 
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        # file_name = "7_3_layer_cropped.vtk"
        # base_name_out = "vessel"

        # process_vtk_file(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

        file_name = "12_three_layer_192.stl"
        base_name_out = "vessel_192"

        process_vtk_file_pyvista(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)


        # must use ring parser file for these 
        # inlet outlet files are vertex files only 
        # spring_strength_rel = 0.0
        # target_strength  = 0.0
        # damping_strength = 0.0
        
        # file_name = "right_ventricle.vtk"
        # base_name_out = "right_ventricle_bdry"
        # process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)            

        # file_name = "right_pa.vtk"
        # base_name_out = "right_pa_bdry"
        # process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)

        # file_name = "left_pa.vtk"
        # base_name_out = "left_pa_bdry"
        # process_vtk_file(file_name, base_name_out, spring_strength_rel, target_strength, damping_strength, scaling)

    mri4d_compare_aortic_files_384 = False 
    if mri4d_compare_aortic_files_384:

        # these are already in CGS units 
        scaling = 1.0

        # target strength aortic_
        target_strength = 97049.24295364547288

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.0147
        kappa_abs = 0.01 * 4.2649e+06 / 3.0 
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        # file_name = "10_3_layer_pt025cm_cropped.vtk"
        # file_name = "8_3_layer_pt052cm.vtk"
        file_name = "13_3_layer_dspt025_five_layer_inlet.vtk"
        base_name_out = "vessel_384"

        process_vtk_file(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

        base_name_out = "vessel_pyvista_384"
        process_vtk_file_pyvista(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

    mri4d_compare_aortic_files_768 = False 
    if mri4d_compare_aortic_files_768:

        # these are already in CGS units 
        scaling = 1.0

        # target strength aortic_192
        target_strength = 388196.97181458189152 * 0.25 * 0.25 

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.014

        kappa_abs = 0.25 * 0.01 * 4.2649e+06 / 3.0 
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        file_name = "12_three_layer_fix_extender_mesh_pt02_ext_pt014.stl"
        base_name_out = "vessel_768"

        process_vtk_file_pyvista(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)


    aorta_192_files = False 
    if aorta_192_files:

        # these are already in CGS units 
        scaling = 1.0

        # target strength aortic_192
        target_strength = 232918.18308874915238

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.05
        kappa_abs = 543733.53989471553359  # abs spring constant from cylinder mesh 
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        file_name = "aorta_192.vtk"
        base_name_out = "aorta_192"

        process_vtk_file(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

    aorta_384_files = False 
    if aorta_384_files:

        # these are already in CGS units 
        scaling = 1.0

        # target strength aortic_384
        target_strength = 58229.54577218728809

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.025

        # abs spring constant from cylinder mesh 
        # scale this down by 2, and down by 2 again with ds
        # if a single spring is split its kappa_abs goes down by two 
        # but there are more of them, if not in a regular way 
        # could turn down by 4 if needed for stability
        kappa_abs = 543733.53989471553359 / 2  
        
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        file_name = "aorta_384.vtk"
        base_name_out = "aorta_384"

        process_vtk_file(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)


    aorta_384_files_bicuspidization = True
    if aorta_384_files_bicuspidization:

        # this is in mm 
        scaling = 0.1

        # target strength aortic_384
        target_strength = 2.0 * 58229.54577218728809

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.025

        # abs spring constant from cylinder mesh 
        # scale this down by 2, and down by 2 again with ds
        # if a single spring is split its kappa_abs goes down by two 
        # but there are more of them, if not in a regular way 
        # could turn down by 4 if needed for stability
        kappa_abs = 0.01 * 543733.53989471553359 / 2  
        
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        # file_name = "4_aorta_remeshed_pt25mm_3cm_extender_layers_constriction.stl"
        file_name = "6_aorta_remeshed_pt25mm_3cm_extenders_layers_double_constriction.stl"
        base_name_out = "aorta_384"

        process_vtk_file_pyvista(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

    aorta_192_files_bicuspidization = False
    if aorta_192_files_bicuspidization:

        # this is in mm 
        scaling = 0.1

        # target strength aortic_192
        target_strength = 232918.18308874915238

        # absolute spring const for cross layer springs of length 
        ds_extrude = 0.05
        kappa_abs = 0.01 * 543733.53989471553359  # abs spring constant from cylinder mesh 
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0

        # file_name = "22_meshed_pt5mm_3_layer_5_layer_inlet_constriction.stl"
        # file_name = "8_5_layer_pt5mm_lv_cropped.stl"
        # file_name = "9_aorta_remeshed_pt5mm_2_cm_extender_5_layer_constriction.stl"
        # file_name = '6_aorta_remeshed_pt5mm_2cm_extender_layers_constriction.stl'
        file_name = '7_aorta_remeshed_pt5mm_2cm_extender_layers_double_constriction.stl'

        base_name_out = "aorta_192"

        process_vtk_file_pyvista(file_name, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)

