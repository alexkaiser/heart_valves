import os 
import pyvista 
import numpy as np 



def extra_radius(x, x_min, x_max, extension_value):

    if x < x_min:
        return 0.0
    elif x < x_max:
        # 0*(x - x_max)/(x_min - x_max) + extension_value*(x - x_min)/(x_max - x_min)
        # linearly interpolate from 0 to extension_value
        return extension_value*(x - x_min)/(x_max - x_min)
    else:
        return extension_value
    

def morph_extender(mesh, 
                   fname_out, 
                   mesh_boundary, 
                   normal_direction, 
                   extension_value,
                   masking_width, 
                   z_translation_max = 0.0,
                   enforce_flat_bdry = True, 
                   flat_bdry_tolerance = 1.0e-3):


    if normal_direction != 0:
        raise NotImplmentedError("normal must be x direction for now")

    # max to mask over 
    x_max = np.max(mesh.points[:,normal_direction])
    x_min = x_max - masking_width
    print("x_min = ", x_min, "x_max = ", x_max)

    # compute the centroid of the mesh 
    centroid = np.mean(mesh_boundary.points, axis=0)
    print("centroid = ", centroid)

    mesh_adjusted = mesh 

    for idx, pt in enumerate(mesh.points):

        normal = (pt - centroid)
        normal[0] = 0.0
        normal /= np.linalg.norm(normal) 

        extra_r = extra_radius(pt[0], x_min, x_max, extension_value)

        z_inc = extra_radius(pt[0], x_min, x_max, z_translation_max)


        increment = [0, extra_r * normal[1], extra_r * normal[2] - z_inc]

        if enforce_flat_bdry:
            if abs(pt[0] - x_max) < flat_bdry_tolerance:
                pt[0] = x_max

        mesh_adjusted.points[idx] = pt + increment 

    mesh_adjusted.save(fname_out)
    return mesh_adjusted


if __name__== "__main__":

    res_192 = True 
    if res_192:
        fname = "21_meshed_pt5mm_3_layer_5_layer_inlet.stl"
        fname_out = "22_meshed_pt5mm_3_layer_5_layer_inlet_constriction.stl"

        bdry_filename = 'lvot_bdry_192.vtu'
        bdry_filename_out = 'lvot_bdry_192_constriction.vtu'

    else:
        # basic 
        fname = "16_meshed_pt25mm_3_layer_5_layer_inlet_outlet.stl"
        fname_out = "17_meshed_pt25mm_3_layer_5_layer_inlet_outlet_constriction.stl"

        bdry_filename = 'lvot_bdry_384_layer_3_pt5mm.vtu'
        bdry_filename_out = 'lvot_bdry_384_layer_3_pt5mm_constriction.vtu'



    
    mesh = pyvista.read(fname)
    mesh_boundary = pyvista.read(bdry_filename)

    # x direction 
    normal_direction = 0

    # mesh in mm, mask 1 cm worth 
    masking_width = 10.0

    # if true, adjusts x component to be exactly equal to this value 
    enforce_flat_bdry = True
    flat_bdry_tolerance = 1.0e-3

    # 1 mm out at the ends 
    extension_value = 0.0

    z_translation_max = 5.0

    mesh_adjusted = morph_extender(mesh, 
                                   fname_out, 
                                   mesh_boundary, 
                                   normal_direction, 
                                   extension_value,
                                   masking_width, 
                                   z_translation_max,
                                   enforce_flat_bdry, 
                                   flat_bdry_tolerance)

    # pyvista.plot(mesh_adjusted)

    mesh_boundary_copy = mesh_boundary

    mesh_boundary_adjusted = morph_extender(mesh_boundary_copy, 
                                            bdry_filename_out, 
                                            mesh_boundary, 
                                            normal_direction, 
                                            extension_value,
                                            masking_width,
                                            z_translation_max, 
                                            enforce_flat_bdry, 
                                            flat_bdry_tolerance)

    # pyvista.plot(mesh_boundary_adjusted)

    



