
import os 
import pyvista 
import numpy as np 
import math 
import pdb


def add_layer(mesh, dist):

    mesh.compute_normals(inplace=True)

    mesh_extruded = mesh.copy()

    mesh_extruded.points += dist * mesh_extruded.point_data['Normals']

    return mesh_extruded


def make_layered_mesh(mesh, dist, n_layers=3, crop_layers=0, crop_dist=1.0, crop_direction=None):

    mesh_layers = mesh.copy()

    for layer in range(1,n_layers):

        dist_temp = dist * layer
        mesh_extruded = add_layer(mesh, dist_temp)

        mesh_layers = pyvista.merge([mesh_layers, mesh_extruded])

    if crop_layers > 0:

        if crop_direction is None:
            raise ValueError('must provide crop_direction if crop_layers > 0')

        if len(crop_direction) != 6:
            raise ValueError('crop_direction must be len 6 vector')

        # format is (xmin, xmax, ymin, ymax, zmin, zmax)
        bounds = np.array(mesh.bounds)

        print("bounds original = ", bounds)

        bounds -= crop_dist * crop_direction

        print("bounds adjusted = ", bounds)

        for layer in range(crop_layers):

            mesh_clipped = mesh.clip_box(bounds).extract_surface()

            print("mesh_clipped = ", mesh_clipped)

            dist_temp = dist * (layer + n_layers)
            mesh_extruded = add_layer(mesh_clipped, dist_temp)            

            mesh_layers = pyvista.merge([mesh_layers, mesh_extruded])


    return mesh_layers

def extract_boundary_meshes(mesh, dist, layer_num, name_largest, name_smallest):

    mesh_layers = mesh.copy()
    dist_temp = dist * layer_num
    mesh_extruded = add_layer(mesh, dist_temp)

    mesh_edges = mesh_extruded.extract_feature_edges(boundary_edges=True,
                                                     non_manifold_edges=False,
                                                     feature_edges=False,
                                                     manifold_edges=False)    

    # quick santity check for number of regions 
    conn = mesh_edges.connectivity()
    assert np.max(conn.point_data['RegionId']) <= 1

    mesh_largest = pyvista.UnstructuredGrid(mesh_edges.connectivity('specified', region_ids=[0]))
    mesh_largest.save(name_largest)

    mesh_smallest = pyvista.UnstructuredGrid(mesh_edges.connectivity('specified', region_ids=[1]))
    mesh_smallest.save(name_smallest)



def mesh_info(mesh, dx, scaling=1.0):

    mesh.points *= scaling

    print("mesh.bounds = ", mesh.bounds)

    NX = (mesh.bounds[1] - mesh.bounds[0])/dx
    NY = (mesh.bounds[3] - mesh.bounds[2])/dx
    NZ = (mesh.bounds[5] - mesh.bounds[4])/dx

    print("NX NY NZ = ", NX, ", ", NY, ", ", NZ)

    ymid = np.mean(mesh.bounds[2:4])

    print("ymid = ", ymid)




if __name__== "__main__":

    filename = '3_aorta_remeshed_pt25_3cm_extenders_no_cap.stl'

    filename_out = '4_aorta_remeshed_pt25_3cm_extenders_3_layer_5_layer_extender.stl'

    mesh = pyvista.read(filename)

    dist = .25

    n_layers = 3

    # mesh_extruded = add_layer(mesh, dist)

    crop_layers = 2
    crop_dist = 30.0
    # crop x,z maximum 
    crop_direction = np.array([0, 1, 0, 0, 0, 1])

    mesh_layers = make_layered_mesh(mesh, dist, n_layers, crop_layers, crop_dist, crop_direction)

    # print summary 
    N = 192
    dx_384 = 0.1 * (192.0/N)
    scaling = 0.1
    mesh_info(mesh_layers, dx_384, scaling)
    quit() 


    mesh_layers.save(filename_out)

    # note zero indexing on layers 
    layer_num_edges = 2

    name_largest = 'lvot_boundary_layer_3.vtu'
    name_smallest = 'aorta_boundary_layer_3.vtu'

    extract_boundary_meshes(mesh, dist, layer_num_edges, name_largest, name_smallest)


