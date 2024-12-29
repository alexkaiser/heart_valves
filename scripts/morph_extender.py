import os 
import pyvista 
import numpy as np 
import math 
import pdb


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
                   fname_out, 
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

    mesh_adjusted.save(fname_out)
    return mesh_adjusted


if __name__== "__main__":

    res_192 = True
    if res_192:
        fname = "5_aorta_remeshed_pt5mm_2_cm_extender_layers.stl"
        # fname_out = "6_aorta_remeshed_pt5mm_2cm_extender_layers_constriction.stl"
        fname_out = "7_aorta_remeshed_pt5mm_2cm_extender_layers_double_constriction.stl"

        bdry_filename = 'lvot_bdry.vtu'
        bdry_filename_out = 'lvot_bdry_192_layer_3_constriction.vtu'

        bdry_filename_aorta = 'aorta_bdry.vtu'
        bdry_filename_aorta_out = 'aorta_bdry_192_layer_3_double_constriction.vtu'

        # bdry_filename = "lvot_bdry_192_layer_1.vtu"
        # bdry_filename_out = 'lvot_bdry_384_layer_1_pt5mm_constriction.vtu'

    else:
        # basic 
        fname = "3_aorta_remeshed_pt25mm_3cm_extender_layers.stl"
        fname_out = "4_aorta_remeshed_pt25mm_3cm_extender_layers_constriction.stl"

        bdry_filename = 'lvot_bdry.vtu'
        bdry_filename_out = 'lvot_bdry_384_layer_3_constriction.vtu'



    
    mesh = pyvista.read(fname)
    mesh_boundary = pyvista.read(bdry_filename)

    mesh_boundary_aorta = pyvista.read(bdry_filename_aorta)

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
                                   fname_out, 
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
                                            bdry_filename_out, 
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
    masking_width = 10.0
    mesh_adjusted = morph_extender(mesh_adjusted, 
                                   fname_out, 
                                   mesh_boundary_aorta, 
                                   normal_direction, 
                                   extension_value,
                                   masking_width, 
                                   enforce_flat_bdry, 
                                   flat_bdry_tolerance, 
                                   cos_interpolation)        

    mesh_aorta_boundary_adjusted = morph_extender(mesh_boundary_aorta, 
                                            bdry_filename_aorta_out, 
                                            mesh_boundary_aorta, 
                                            normal_direction, 
                                            extension_value,
                                            masking_width,
                                            enforce_flat_bdry, 
                                            flat_bdry_tolerance,
                                            cos_interpolation)    



