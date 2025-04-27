import pyvista 
import os 
import numpy as np 


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



if __name__== "__main__":
                
    run_192 = False
    run_384_shell = True
    if run_192:
        fname_in = "4_aorta_remeshed_pt5mm_2_cm_extender.stl"
        fname_out = "5_aorta_remeshed_pt5mm_2_cm_extender_layers.stl"

        n_layers_full = 3
        n_layers_extenders = 2

        # extrude length in mm 
        ds = 0.5

    elif run_384_shell:
        fname_in = "3_aorta_remeshed_pt25_3cm_extenders_no_cap.stl"
        fname_out = "3_aorta_remeshed_pt25_3cm_extenders_no_cap_two_layer.stl"

        n_layers_full = 2
        n_layers_extenders = 0

        # extrude length in mm 
        ds = 0.5

    else:
        fname_in = "2_aorta_remeshed_pt25mm_3cm_extender.stl"
        fname_out = "3_aorta_remeshed_pt25mm_3cm_extender_layers.stl"

        n_layers_full = 3
        n_layers_extenders = 2

        # extrude length in mm 
        ds = 0.5/2

    mesh = pyvista.read(fname_in)


    extender_direction_idx = [0,2]
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

    mesh_combined.save(fname_out)

    boundary_meshes = False

    if boundary_meshes:

        z_max = np.max(mesh_combined.points[:,2])
        tol_edges = 1.0e-1

        print(f"x_max = {x_max:.16f}")
        print(f"z_max = {z_max:.16f}")

        # write out edge meshes 
        edge_name_const_x = "lvot_bdry.vtu"
        edge_name_const_z = "aorta_bdry.vtu"

        conn = edges.connectivity()

        print("conn = ", conn, "conn.point_data = ", conn.point_data)

        # n_regions = 2 

        n_regions = max(conn.point_data['RegionId'] + 1)

        min_points_valid_bdry = 50

        # if n_regions != 2:
        #     print("n_regions =", n_regions)
        #     print("conn.point_data['RegionId'] = ", conn.point_data['RegionId'])
        #     raise ValueError("must have two regions")

        for region_id in range(n_regions):

            name_found = False

            bdry_mesh = conn.threshold([region_id, region_id + 0.5], scalars="RegionId", all_scalars=True)

            # print("bdry_mesh = ", bdry_mesh, "bdry_mesh.points = ", bdry_mesh.points)
            print("bdry_mesh = ", bdry_mesh)

            if bdry_mesh.n_points > min_points_valid_bdry:

                if all(abs(bdry_mesh.points[:,0] - x_max) < tol_edges):
                    bdry_mesh_name = edge_name_const_x
                    name_found = True 
                elif all(abs(bdry_mesh.points[:,2] - z_max) < tol_edges):
                    bdry_mesh_name = edge_name_const_z
                    name_found = True 

                if name_found:
                    bdry_mesh.save(bdry_mesh_name)
                else:
                    print("Did not find boundary mesh for RegionId ", region_id, ", skipping.")

            else: 
                print("RegionId ", region_id, " is small and probably an artifact with ", bdry_mesh.n_points, " points, skipping.")



