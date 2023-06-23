import sys
import os 
import shutil 
import pyvista
import numpy as np 
import glob


def convert_csv(basename, suffix='_faces', vertex_ext='_vertices.csv', cells_ext='_cells.csv', ext_out='.stl'):

    vertices = np.loadtxt(basename + vertex_ext)

    cells = np.loadtxt(basename + cells_ext, dtype=int)

    # print("vertices = ", vertices)
    # print("cells = ", cells)

    n_elem = cells.shape[0]

    # cell_type = np.empty(n_elem, np.uint8)
    cell_type = np.empty(n_elem, np.uint)
    cell_type[:] = [pyvista.CellType.HEXAHEDRON]

    mesh = pyvista.UnstructuredGrid(cells, cell_type, vertices)

    # mesh.plot(show_edges=True)

    outname = basename + suffix + ext_out

    print("saving outname = ", outname)

    if (ext_out == '.stl') or (ext_out == '.vtp'):
        mesh_surface = mesh.extract_surface()
        mesh_surface.save(outname)
    else:
        # volumetric default 
        mesh.save(outname)





if __name__ == '__main__':


    run_all = True 

    if run_all:

        basename_no_frame = 'aortic_no_partition_384'
        suffix='_faces'
        vertex_ext = '_vertices.csv'
        extension_out = '.vtu'

        file_list = glob.glob(basename_no_frame + '*' + vertex_ext)

        for name_full in file_list:

            basename = name_full.rsplit(vertex_ext)[0]
            convert_csv(basename, suffix, vertex_ext=vertex_ext, ext_out=extension_out)


    local_example = False 

    if local_example:

        # dir_list = ["aortic_14234920_384_a75f53c_true_bicuspid_circ_2pt8_1pt1d_rad_1pt4_new_initial_cond", 
        #     "aortic_14103984_384_28eb02c_true_bicuspid_circ_3pt0_1pt2d_rad_1pt7_new_initial_cond", 
        #     "aortic_14097127_384_4b0882c_true_bicuspid_circ_3pt5_1pt4d_rad_1pt7_new_initial_cond_2", 
        #     "aortic_14095958_384_67bdcde_true_bicuspid_circ_3pt9_1pt65d_rad_1pt7_new_initial_cond_2", 
        #     "aortic_14149904_384_a75f53c_true_bicuspid_circ_4pt5_1pt8d_rad_1pt7_new_initial_cond", 
        #     "aortic_14728686_384_a75f53c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt7_circ_prestretch_1pt08", 
        #     "aortic_15051277_384_98d4f6c_true_bicuspid_circ_3pt9_1pt57d_rad_1pt4", 
        #     "aortic_16698035_384_23a9213_true_bicuspid_circ_3pt9_1pt6d_rad_1pt65_less_bowl"]; 

        # suffixes = ["_circ_1pt1d_rad_1pt4", 
        #     "_circ_1pt2d_rad_1pt7", 
        #     "_circ_1pt4d_rad_1pt7", 
        #     "_circ_1pt57d_rad_1pt7", 
        #     "_circ_1pt8d_rad_1pt7",
        #     "_circ_1pt57d_rad_1pt7_circ_prestretch_1pt08", 
        #     "_circ_1pt57d_rad_1pt4", 
        #     "_circ_1pt6d_rad_1pt65_less_bowl"]; 

        suffixes = ['circ_3pt0_1pt2d_rad_1pt4_less_bowl']; 

        dir_list = ['aortic_21254567_384_fc9c4b3_true_bicuspid_circ_3pt0_1pt2d_rad_1pt4_less_bowl_layers_2e4']; 


        base_dir = os.getcwd()

        for direc, suffix in zip(dir_list, suffixes):

            print("processing direc ", direc)

            os.chdir(direc)

            print("in         direc ", os.getcwd())

            frame_numbers = [241, 406]

            # frame_numbers = [241]

            basename_no_frame = 'aortic_no_partition_384'

            extension_out = '.vtu'

            for frame in frame_numbers:

                basename = basename_no_frame + str(frame).zfill(4)

                convert_csv(basename, suffix, ext_out=extension_out)

            os.chdir(base_dir)







