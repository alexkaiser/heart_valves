import pyvista 
import os 
import numpy as np 
import math 

def expand_mesh(mesh, 
                ds,
                n_layers_full=3, 
                n_layers_extenders=0,
                extender_direction_idx=None, 
                extender_top=True,
                extender_width=1.0,
                extract_edge_layer=None):

    mesh.compute_normals(cell_normals=False, point_normals=True, split_vertices=True, inplace=True, feature_angle=70)

    mesh_combined = mesh

    edges = None  

    for i in range(1,n_layers_full):

        # copy original layer 
        mesh_new_layer = mesh.copy(deep=True)

        mesh_new_layer.points += ds * i * mesh.point_data['Normals']

        # print("mesh_new_layer = ", mesh_new_layer)

        mesh_combined = mesh_combined.merge(mesh_new_layer, merge_points=False)
        
        if i == extract_edge_layer:
            edges = mesh_new_layer.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)

        # print("mesh_combined = ", mesh_combined)


    if n_layers_extenders > 0:

        if extender_direction_idx is None: 
            raise ValueError("must provide extender_directions if requesting extenders")

        if not extender_top:
            raise NotImplementedError("top extenders only implemented")

        for direction, width in zip(extender_direction_idx, extender_width):

            if direction not in range(3):
                raise ValueError()

            for i in range(n_layers_extenders):
                # extrude whole mesh 
                mesh_new_layer = mesh.copy(deep=True)
                mesh_new_layer.points += ds * (i + n_layers_full) * mesh.point_data['Normals']
            
                # clip the layer 
                max_val = np.max(mesh.points[:,direction])
                crop_location = max_val - width

                direction_vector = tuple(-1.0 * (j == direction) for j in range(3))
                origin = tuple(crop_location * (j == direction) for j in range(3))

                mesh_new_layer.clip(normal=direction_vector, 
                                    origin=origin, 
                                    inplace=True)

                mesh_combined = mesh_combined.merge(mesh_new_layer, merge_points=False)


    return mesh_combined, edges


def find_boundary_meshes(mesh, 
                         edges,
                         tol_edges_rel = 1.0e-3,
                         min_points_valid_bdry = 50,
                         debug_output = True, 
                         enforce_flat = True):
    '''
    Find all boundary meshes among connected component of edge meshs 
    '''

    boundary_meshes = []
    components = []
    sides = []

    conn = edges.connectivity()

    if debug_output:
        print("conn = ", conn, "conn.point_data = ", conn.point_data)

    n_regions = max(conn.point_data['RegionId'] + 1)
    print("found n_regions = ", n_regions)

    for region_id in range(n_regions):
        print("testing region ", region_id)

        bdry_mesh = conn.threshold([region_id, region_id + 0.5], scalars="RegionId", all_scalars=True)

        print("bdry_mesh.n_points = ", bdry_mesh.n_points)

        if bdry_mesh.n_points > min_points_valid_bdry:

            for component in range(3):

                box_width = mesh_combined.bounds[2*component + 1] - mesh_combined.bounds[2*component]

                for side in range(2):

                    # mesh_combined.bounds (x_min, x_max, y_min, y_max, z_min, z_max)
                    boundary_value = mesh_combined.bounds[2*component + side]

                    if all(abs(bdry_mesh.points[:,component] - boundary_value)/box_width < tol_edges_rel):

                        if debug_output:
                            print("Found mesh on component ", component, ", side ", side)

                        if enforce_flat:
                            for i in range(bdry_mesh.n_points):
                                bdry_mesh.points[i,component] = boundary_value

                            if max(abs(bdry_mesh.points[:,component] - boundary_value) > 0):
                                print("max diff = ", max(abs(bdry_mesh.points[:,component] - boundary_value)))
                                raise ValueError("Flat points enforced but not actually flat")

                        boundary_meshes.append(bdry_mesh)

        else: 
            if debug_output:
                print("RegionId ", region_id, " is small and probably an artifact with ", bdry_mesh.n_points, " points, skipping.")

    return boundary_meshes


def extra_radius(x, x_min, x_max, extension_value, cos_interpolation=False):

    if x < x_min:
        return 0.0
    elif x < x_max:
        # 0*(x - x_max)/(x_min - x_max) + extension_value*(x - x_min)/(x_max - x_min)
        # linearly interpolate from 0 to extension_value
        if cos_interpolation:
            arg = (x - x_min)/(x_max - x_min)
            return extension_value * 0.5*(-math.cos(math.pi * arg) + 1.0)
        else:   
            return extension_value*(x - x_min)/(x_max - x_min)
    else:
        return extension_value
    

def morph_extender(mesh, 
                   mesh_boundary, 
                   normal_direction, 
                   extension_value,
                   masking_width, 
                   enforce_flat_bdry = True, 
                   flat_bdry_tolerance = 1.0e-3,
                   cos_interpolation = False):


    # max to mask over 
    pt_normal_max = np.max(mesh.points[:,normal_direction])
    pt_normal_min = pt_normal_max - masking_width
    print("pt_normal_min = ", pt_normal_min, "pt_normal_max = ", pt_normal_max)

    # compute the centroid of the mesh 
    centroid = np.mean(mesh_boundary.points, axis=0)
    print("centroid = ", centroid)

    mesh_adjusted = mesh 

    for idx, pt in enumerate(mesh.points):

        normal = (pt - centroid)
        normal[normal_direction] = 0.0
        normal /= np.linalg.norm(normal) 

        extra_r = extra_radius(pt[normal_direction], pt_normal_min, pt_normal_max, extension_value, cos_interpolation)

        increment = extra_r * normal

        if enforce_flat_bdry:
            if abs(pt[normal_direction] - pt_normal_max) < flat_bdry_tolerance:
                pt[normal_direction] = pt_normal_max

        mesh_adjusted.points[idx] = pt + increment 

    return mesh_adjusted


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




def process_mesh_to_vertex_spring_target(mesh, base_name_out, spring_strength_rel, target_strength, damping_strength=0.0, scaling=1.0, zero_springs=False):

    print("meah.n_points = ", mesh.n_points)

    mesh_unst = pyvista.UnstructuredGrid(mesh)
    # mesh_unst.save(base_name_out + ".vtu")

    print("mesh_unst.n_points = ", mesh_unst.n_points)

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


def process_ring_to_vertex(mesh, base_name_out, scaling=1.0):

    vertex_file = open(base_name_out + '.vertex', 'w')

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

def print_boundary_info(mesh, DX, extra_bdry_min=0.25, multiple=8):

    x_min = mesh.bounds[0]
    x_max = mesh.bounds[1]
    y_min = mesh.bounds[2]
    y_max = mesh.bounds[3]
    z_min = mesh.bounds[4]
    z_max = mesh.bounds[5]

    NX_temp = (x_max - x_min)/DX
    NY_temp = (y_max - y_min)/DX
    NZ_temp = (z_max - z_min)/DX

    # ensure all are muliples of prescribed number for multigrid
    NX = multiple * math.ceil(NX_temp/multiple)
    NY = multiple * math.ceil(NY_temp/multiple)
    NZ = multiple * math.ceil(NZ_temp/multiple)

    X_HIGH = x_max
    Y_MID  = 0.5 * (y_max + y_min)
    Y_HIGH = Y_MID + (NY*DX)/2.0
    Z_HIGH = z_max

    X_LOW = X_HIGH - DX*NX
    Y_LOW = Y_HIGH - DX*NY
    Z_LOW = Z_HIGH - DX*NZ

    print("bounds = ", mesh.bounds)
    print("computed boundary = ", X_LOW, " ", X_HIGH, " ", Y_LOW, " ", Y_HIGH, " ", Z_LOW, " ", Z_HIGH)
    print("NX NY NZ = ", NX, " ", NY, " ", NZ)

    print("Before:")
    print("x extra = ", x_min - X_LOW)
    print("y extra bottom = ", y_min - Y_LOW)
    print("y extra top    = ", Y_HIGH - y_max)
    print("z extra        = ", z_min - Z_LOW)

    if (x_min - X_LOW) < extra_bdry_min:
        NX += multiple
        X_LOW = X_HIGH - DX*NX
        if (x_min - X_LOW) < extra_bdry_min:
            raise ValueError('did not fix minimum boundary')

    if (y_min - Y_LOW) < extra_bdry_min:
        NY += multiple
        Y_HIGH = Y_MID + (NY*DX)/2.0
        Y_LOW = Y_HIGH - DX*NY
        if (y_min - Y_LOW) < extra_bdry_min:
            raise ValueError('did not fix minimum boundary')

    if (Y_HIGH - y_max) < extra_bdry_min:
        NY += multiple
        Y_HIGH = Y_MID + (NY*DX)/2.0
        Y_LOW = Y_HIGH - DX*NY
        if (Y_HIGH - y_max) < extra_bdry_min:
            raise ValueError('did not fix maximum boundary')

    if (z_min - Z_LOW) < extra_bdry_min:
        NZ += multiple
        Z_LOW = Z_HIGH - DX*NZ
        if (z_min - Z_LOW) < extra_bdry_min:
            raise ValueError('did not fix minimum boundary')
    
    print("After:")
    print("x extra = ", x_min - X_LOW)
    print("y extra bottom = ", y_min - Y_LOW)
    print("y extra top    = ", Y_HIGH - y_max)
    print("z extra        = ", z_min - Z_LOW)

    print("computed boundary = ", X_LOW, " ", X_HIGH, " ", Y_LOW, " ", Y_HIGH, " ", Z_LOW, " ", Z_HIGH)
    print("NX NY NZ = ", NX, " ", NY, " ", NZ)

    

    print(f"X_HIGH = {X_HIGH:.14f}")
    print(f"Z_HIGH = {Z_HIGH:.14f}")
    print(f"Y_MID  = {Y_MID:.14f}")

    print("NX = ", NX)
    print("NY = ", NY)
    print("NZ = ", NZ)


if __name__== "__main__":
    
    example = False 
    if example:

        fname_in = "2_aorta_lv_extender.stl"
        fname_out = "3_aorta_lv_extender_layers.stl"

        n_layers_full = 3
        n_layers_extenders = 2

        # extrude length in mm 
        ds = 0.5        

        mesh = pyvista.read(fname_in)


        extender_direction_idx = [0,2] # extra mesh layers at inlet and outlet 
        extender_top = True
        extender_width = [30.0, 10.0]
        extract_edge_layer = 2

        mesh_combined, edges = expand_mesh(mesh,
                                           ds, 
                                           n_layers_full, 
                                           n_layers_extenders, 
                                           extender_direction_idx,
                                           extender_top,
                                           extender_width,
                                           extract_edge_layer)

        boundary_meshes = find_boundary_meshes(mesh_combined, edges)

        # for idx, bdry_mesh in enumerate(boundary_meshes):
        #     bdry_mesh.save('boundary_mesh' + str(idx).zfill(4) + '.vtu')

        # 
        mesh = mesh_combined
        mesh_boundary = boundary_meshes[0]
        mesh_boundary_aorta = boundary_meshes[1]

        # x direction 
        normal_direction = 0

        # mesh in mm
        masking_width = 15.0

        # if true, adjusts x component to be exactly equal to this value 
        enforce_flat_bdry = True
        flat_bdry_tolerance = 1.0e-3

        # 1 mm out at the ends 
        extension_value = 5.0

        cos_interpolation = True

        mesh_adjusted = morph_extender(mesh, 
                                       mesh_boundary, 
                                       normal_direction, 
                                       extension_value,
                                       masking_width, 
                                       enforce_flat_bdry, 
                                       flat_bdry_tolerance, 
                                       cos_interpolation)

        # pyvista.plot(mesh_adjusted)

        mesh_boundary_copy = mesh_boundary

        mesh_boundary_adjusted = morph_extender(mesh_boundary_copy, 
                                                mesh_boundary, 
                                                normal_direction, 
                                                extension_value,
                                                masking_width,
                                                enforce_flat_bdry, 
                                                flat_bdry_tolerance,
                                                cos_interpolation)

        # pyvista.plot(mesh_boundary_adjusted)

        # aorta side 
        normal_direction = 2
        masking_width = 15.0
        mesh_adjusted = morph_extender(mesh_adjusted, 
                                       mesh_boundary_aorta, 
                                       normal_direction, 
                                       extension_value,
                                       masking_width, 
                                       enforce_flat_bdry, 
                                       flat_bdry_tolerance, 
                                       cos_interpolation)        

        mesh_aorta_boundary_adjusted = morph_extender(mesh_boundary_aorta, 
                                                mesh_boundary_aorta, 
                                                normal_direction, 
                                                extension_value,
                                                masking_width,
                                                enforce_flat_bdry, 
                                                flat_bdry_tolerance,
                                                cos_interpolation)    


        mesh_adjusted.save("vessel_post_morph.stl")
        mesh_boundary_adjusted.save("lvot_bdry_morph.vtu")
        mesh_aorta_boundary_adjusted.save("aorta_bdry_morph.vtu")


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

        base_name_out = "aorta_384"

        process_mesh_to_vertex_spring_target(mesh_adjusted, base_name_out, kappa_rel, target_strength, damping_strength, scaling, zero_springs)


        scaling = 0.1

        base_name_out = "aorta_bdry_384"
        process_ring_to_vertex(mesh_aorta_boundary_adjusted, base_name_out, scaling)

        base_name_out = "lvot_bdry_384"
        process_ring_to_vertex(mesh_boundary_adjusted, base_name_out, scaling)


    historical_3 = True 
    if historical_3:

        fine_res = True
        if fine_res:

            native = False 
            z_zero = False 
            z_minus_two = True 

            if native:
                fname_in = "2_aorta_native_original_cropped_pt25.stl"
                fname_out = "aorta_native_hist3_pt25mm_384.stl"
                fname_out_unstructured = "aorta_native_hist3_pt25mm_384.vtu"
                aorta_name = "aorta_native_hist3_pt25mm_384"
                ao_boundary_name = "aorta_native_hist3_bdry_384"
                lvot_boundary_name = "lvot_native_hist3_bdry_384"
            elif z_zero:
                fname_in = "6_aorta_1pt98_original_cropped_pt25.stl"
                fname_out = "aorta_1pt98_hist3_pt25mm_384.stl"
                fname_out_unstructured = "aorta_1pt98_hist3_pt25mm_384.vtu"
                aorta_name = "aorta_1pt98_hist3_pt25mm_384"
                ao_boundary_name = "aorta_1pt98_hist3_bdry_384"
                lvot_boundary_name = "lvot_1pt98_hist3_bdry_384"
            elif z_minus_two:
                fname_in = "8_aorta_1pt7_original_cropped_pt25.stl"
                fname_out = "aorta_1pt7_hist3_pt25mm_384.stl"
                fname_out_unstructured = "aorta_1pt7_hist3_pt25mm_384.vtu"
                aorta_name = "aorta_1pt7_hist3_pt25mm_384"
                ao_boundary_name = "aorta_1pt7_hist3_bdry_384"
                lvot_boundary_name = "lvot_1pt7_hist3_bdry_384"

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

            # extrude length in mm 
            ds = 0.25        

        else:
            fname_in = "4_aorta_native_original_cropped_pt5.stl"
            fname_out = "aorta_native_hist3_pt5mm.stl"
            aorta_name = "aorta_native_hist3_pt5mm_192"
            ao_boundary_name = "aorta_native_hist3_bdry_192"
            lvot_boundary_name = "lvot_native_hist3_bdry_192"

            # target strength aortic_192
            target_strength = 232918.18308874915238

            # absolute spring const for cross layer springs of length 
            ds_extrude = 0.05
            kappa_abs = 0.01 * 543733.53989471553359  # abs spring constant from cylinder mesh 

            # extrude length in mm 
            ds = 0.5        

        n_layers_full = 3
        n_layers_extenders = 2

        scaling = 0.1


        mesh = pyvista.read(fname_in)

        extender_direction_idx = [0,2] # extra mesh layers at inlet and outlet 
        extender_top = True
        extender_width = [30.0, 30.0]
        extract_edge_layer = 2

        mesh_combined, edges = expand_mesh(mesh,
                                           ds, 
                                           n_layers_full, 
                                           n_layers_extenders, 
                                           extender_direction_idx,
                                           extender_top,
                                           extender_width,
                                           extract_edge_layer)

        print("mesh_combined.n_points = ", mesh_combined.n_points)

        boundary_meshes = find_boundary_meshes(mesh_combined, edges)

        print("boundary_meshes = ", boundary_meshes)

        # for idx, bdry_mesh in enumerate(boundary_meshes):
        #     bdry_mesh.save('boundary_mesh' + str(idx).zfill(4) + '.vtu')

        # 
        mesh = mesh_combined
        mesh_boundary = boundary_meshes[0]
        mesh_boundary_aorta = boundary_meshes[1]

        # x direction 
        normal_direction = 0

        # mesh in mm
        masking_width = 15.0

        # if true, adjusts x component to be exactly equal to this value 
        enforce_flat_bdry = True
        flat_bdry_tolerance = 1.0e-3

        # 1 mm out at the ends 
        extension_value = 5.0

        cos_interpolation = True

        mesh_adjusted = morph_extender(mesh, 
                                       mesh_boundary, 
                                       normal_direction, 
                                       extension_value,
                                       masking_width, 
                                       enforce_flat_bdry, 
                                       flat_bdry_tolerance, 
                                       cos_interpolation)

        print("mesh_adjusted.n_points = ", mesh_adjusted.n_points)

        # pyvista.plot(mesh_adjusted)

        mesh_boundary_copy = mesh_boundary

        mesh_boundary_adjusted = morph_extender(mesh_boundary_copy, 
                                                mesh_boundary, 
                                                normal_direction, 
                                                extension_value,
                                                masking_width,
                                                enforce_flat_bdry, 
                                                flat_bdry_tolerance,
                                                cos_interpolation)

        # pyvista.plot(mesh_boundary_adjusted)

        # aorta side 
        normal_direction = 2
        mesh_adjusted = morph_extender(mesh_adjusted, 
                                       mesh_boundary_aorta, 
                                       normal_direction, 
                                       extension_value,
                                       masking_width, 
                                       enforce_flat_bdry, 
                                       flat_bdry_tolerance, 
                                       cos_interpolation)        

        mesh_aorta_boundary_adjusted = morph_extender(mesh_boundary_aorta, 
                                                mesh_boundary_aorta, 
                                                normal_direction, 
                                                extension_value,
                                                masking_width,
                                                enforce_flat_bdry, 
                                                flat_bdry_tolerance,
                                                cos_interpolation)    

        mesh_adjusted.points *= scaling
        mesh_boundary_adjusted.points *= scaling
        mesh_aorta_boundary_adjusted.points *= scaling

        print("before stl save mesh_adjusted.n_points = ", mesh_adjusted.n_points)
        mesh_adjusted.save(fname_out)

        mesh_adjusted_unst = pyvista.UnstructuredGrid(mesh_adjusted)
        print("mesh_adjusted_unst.n_points = ", mesh_adjusted_unst.n_points)
        mesh_adjusted_unst.save(fname_out_unstructured)

        mesh_boundary_adjusted.save(lvot_boundary_name + ".vtu")
        mesh_aorta_boundary_adjusted.save(ao_boundary_name + ".vtu")

        
        kappa_rel = kappa_abs * ds_extrude
        zero_springs = True

        # damping off 
        damping_strength = 0.0


        process_mesh_to_vertex_spring_target(mesh_adjusted, aorta_name, kappa_rel, target_strength, damping_strength)
        process_ring_to_vertex(mesh_aorta_boundary_adjusted, ao_boundary_name)
        process_ring_to_vertex(mesh_boundary_adjusted, lvot_boundary_name)


        dx = 2 * ds * scaling 
        print_boundary_info(mesh_adjusted, dx, extra_bdry_min=0.25, multiple=16)


    test_bounds = True
    if test_bounds:
        mesh = pyvista.read("aorta_native_hist3_pt25mm_384.stl")

        ds = 0.025        
        dx = 2*ds 
        print_boundary_info(mesh, dx, extra_bdry_min=dx*2.5, multiple=8)





        



