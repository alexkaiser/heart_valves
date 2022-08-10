import pyvista 
import os, sys, glob, shutil
import multiprocessing
import pdb
import time 
import math 
import numpy as np
from pykdtree.kdtree import KDTree

def write_pvd(base_name, dt, nsteps, extension, nprocs_sim=1, time_start_1=False):

    prefix = '''<?xml version="1.0"?>
    <VTKFile type="Collection" version="0.1"
             byte_order="LittleEndian"
             compressor="vtkZLibDataCompressor">
      <Collection>
    '''

    suffix = '''  </Collection>
    </VTKFile>
    '''

    initialized = False

    for n in range(nsteps):
        for proc in range(nprocs_sim):

            if not initialized:

                filename_out = base_name + '.pvd'
                print("filename_out = ", filename_out)

                f_write = open(filename_out, 'w')
                f_write.write(prefix)
                initialized = True

            tmp_str = '    <DataSet timestep="'
            if time_start_1:
                tmp_str += '{:.14f}'.format(dt * (n+1))
            else:
                tmp_str += '{:.14f}'.format(dt * n)
            tmp_str += '" group="" part="'
            tmp_str += str(proc) + '"'

            tmp_str += ' file="'
            if nprocs_sim > 1:
                tmp_str += base_name + str(n).zfill(4) + '/' # sorted into directories 
            tmp_str += base_name + str(n).zfill(4) + '.' 
            if nprocs_sim > 1:
                tmp_str += str(proc) + '.'        
            tmp_str += extension
            tmp_str += '"/>\n'

            f_write.write(tmp_str)

    f_write.write(suffix)
    f_write.close()




def interpolate_eulerian_to_mri_mesh(frame_number,
                                     eulerian_basename_orig="eulerian_vars_averaged",
                                     digits_output_eulerian = 4,
                                     eulerian_extension = ".vtu",
                                     suffix_eulerian = '_resampled',
                                     mri_basename = "HealthyNative_",
                                     digits_output_mri = 2,
                                     mri_extension = ".vtk",
                                     output_extension = ".vtu",
                                     boundary_mesh_name = None):

    print("processing frame_number, ", frame_number)

    eulerian_name = eulerian_basename_orig + str(frame_number).zfill(digits_output_eulerian) + eulerian_extension
    eulerian_mesh_orig = pyvista.read(eulerian_name)
    eulerian_mesh_unstructured = pyvista.UnstructuredGrid(eulerian_mesh_orig)
    eulerian_mesh_unstructured_points = eulerian_mesh_unstructured.cell_data_to_point_data()

    mri_name = mri_basename + str(frame_number).zfill(digits_output_mri) + mri_extension
    mri_mesh = pyvista.read(mri_name)
    mri_mesh_unstructured = pyvista.UnstructuredGrid(mri_mesh)

    # restrict MRI mesh to internal 
    if boundary_mesh_name is not None: 

        boundary_mesh = pyvista.read(boundary_mesh_name)
        selected = mri_mesh_unstructured.select_enclosed_points(boundary_mesh, tolerance=0.0, inside_out=False, check_surface=True)

        # all_scalars = True means all points must be within for a cell to stay 
        # all_scalars = False means any point must be within for a cell to stay 
        mri_mesh_unstructured = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

    eulerian_resampled = mri_mesh_unstructured.sample(eulerian_mesh_unstructured_points)


    eulerian_name_out = eulerian_basename_orig + suffix_eulerian
    eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + output_extension

    print("writing file ", eulerian_name_out)
    eulerian_resampled.save(eulerian_name_out)




def run_main():

    run_all = False 

    if run_all:

        # run restricted and not restricted 
        boundary_mesh_names = [None, '2_1_HealthyGeometry_PlanarInletOutlet_smoothed_capped.stl']

        for boundary_mesh_name in boundary_mesh_names:

            nframes = 20

            eulerian_basename_orig="eulerian_vars_averaged"
            digits_output_eulerian = 4
            eulerian_extension = ".vtu"
            suffix_eulerian = '_resampled_points'

            if boundary_mesh_name is not None: 
                suffix_eulerian += '_restricted'
                print("boundary_mesh_name is not None passed")

                # grab this file if it's not here... 
                if not os.path.isfile(boundary_mesh_name):
                    if os.path.isfile(os.path.expanduser('~') + '/mitral_fully_discrete/' + boundary_mesh_name):
                        shutil.copy(os.path.expanduser('~') + '/mitral_fully_discrete/' + boundary_mesh_name, '.') 
                    else: 
                        raise FileNotFoundError("cannot find boundary_mesh_name file = ", boundary_mesh_name)


            print("boundary_mesh_name = ", boundary_mesh_name)
            print("suffix_eulerian = ", suffix_eulerian)

            mri_basename = "HealthyNative_ForSim_"
            digits_output_mri = 2
            mri_extension = ".vtk"
            output_extension = ".vtu"
            
            pool = multiprocessing.Pool(int(nframes/2)) #use all available cores, otherwise specify the number you want as an argument
            for i in range(nframes):
                pool.apply_async(interpolate_eulerian_to_mri_mesh, args=(i,    
                                                                         eulerian_basename_orig,
                                                                         digits_output_eulerian,
                                                                         eulerian_extension,
                                                                         suffix_eulerian,
                                                                         mri_basename,
                                                                         digits_output_mri,
                                                                         mri_extension,
                                                                         output_extension,
                                                                         boundary_mesh_name))
            pool.close()
            pool.join()


            cycle_duration = 8.3250000000000002e-01
            mri_read_times_per_cycle = 10 
            dt_mri_read = cycle_duration / mri_read_times_per_cycle
            output_times_per_cycle   = 20 
            dt_output = cycle_duration / output_times_per_cycle

            time_start_1 = True 

            output_extension = "vtu"
            write_pvd(eulerian_basename_orig + suffix_eulerian, dt_output, output_times_per_cycle, output_extension, time_start_1=time_start_1)



    tests = False
    if tests:

        eulerian_basename_orig = "eulerian_vars_averaged" 
        digits_output_eulerian = 4
        eulerian_extension = ".vtu"
        mri_basename = "HealthyNative_ForSim_"
        digits_output_mri = 2
        mri_extension = ".vtk"

        frame_number = 6 

        eulerian_name = eulerian_basename_orig + str(frame_number).zfill(digits_output_eulerian) + eulerian_extension
        eulerian_mesh_orig = pyvista.read(eulerian_name)
        eulerian_mesh_unstructured = pyvista.UnstructuredGrid(eulerian_mesh_orig)
        eulerian_mesh_unstructured_points = eulerian_mesh_unstructured.cell_data_to_point_data()

        mri_name = mri_basename + str(frame_number).zfill(digits_output_mri) + mri_extension
        mri_mesh = pyvista.read(mri_name)
        mri_mesh_unstructured = pyvista.UnstructuredGrid(mri_mesh)
        # mri_mesh_unstructured_cells = mri_mesh_unstructured.point_data_to_cell_data()

        # not of these working well 
        # eulerian_resampled = mri_mesh_unstructured.sample(eulerian_mesh_unstructured)
        # eulerian_resampled_cells = eulerian_resampled.point_data_to_cell_data()
        # eulerian_resampled_from_mri_cells = mri_mesh_unstructured_cells.sample(eulerian_mesh_unstructured)        
        # eulerian_resampled_from_mri_cells_to_points = eulerian_resampled_from_mri_cells.cell_data_to_point_data()
        
        # onto entire unstructured mesh 
        # eulerian_points_resampled = mri_mesh_unstructured.sample(eulerian_mesh_unstructured_points)

        boundary_mesh_name = '2_1_HealthyGeometry_PlanarInletOutlet_smoothed_capped.stl'
        boundary_mesh = pyvista.read(boundary_mesh_name)
        
        # casting not allowed, need a surfbace mesh 
        # boundary_mesh = pyvista.UnstructuredGrid(boundary_mesh)

        selected = mri_mesh_unstructured.select_enclosed_points(boundary_mesh, tolerance=0.0, inside_out=False, check_surface=True)

        # all_scalars = True means all points must be within for a cell to stay 
        # all_scalars = False means any point must be within for a cell to stay 
        mri_mesh_restricted = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

        eulerian_points_resampled = mri_mesh_restricted.sample(eulerian_mesh_unstructured_points)



        # suffix_eulerian = '_resampled'
        # eulerian_name_out_extension = '.vtu'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_resampled.save(eulerian_name_out)

        # suffix_eulerian = '_resampled_cells'
        # eulerian_name_out_extension = '.vtu'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_resampled_cells.save(eulerian_name_out)

        # suffix_eulerian = '_resampled_from_mri_cells'
        # eulerian_name_out_extension = '.vtu'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_resampled_from_mri_cells.save(eulerian_name_out)

        # suffix_eulerian = '_resampled_from_mri_cells_to_points'
        # eulerian_name_out_extension = '.vtu'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_resampled_from_mri_cells_to_points.save(eulerian_name_out)

        suffix_eulerian = '_points_resampled'
        eulerian_name_out_extension = '.vtu'
        eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        eulerian_points_resampled.save(eulerian_name_out)


        # gaussian at scan radius 
        # does not work written like this, unclear why 
        # radius = 0.09 * 4
        # eulerian_mesh_guass_interp = mri_mesh.interpolate(eulerian_mesh_orig, radius=radius)
        # suffix_eulerian = '_gaussian_interp'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_mesh_guass_interp.save(eulerian_name_out)

    tests_spatial_averaging = True 
    if tests_spatial_averaging:

        eulerian_basename_orig = "eulerian_vars_averaged" 
        digits_output_eulerian = 4
        eulerian_extension = ".vtu"
        mri_basename = "HealthyNative_ForSim_"
        digits_output_mri = 2
        mri_extension = ".vtk"

        frame_number = 6 

        eulerian_name = eulerian_basename_orig + str(frame_number).zfill(digits_output_eulerian) + eulerian_extension
        eulerian_mesh_orig = pyvista.read(eulerian_name)
        eulerian_mesh_unstructured = pyvista.UnstructuredGrid(eulerian_mesh_orig)
        eulerian_mesh_centers = eulerian_mesh_unstructured.cell_centers()

        mri_name = mri_basename + str(frame_number).zfill(digits_output_mri) + mri_extension
        mri_mesh = pyvista.read(mri_name)
        mri_mesh_unstructured = pyvista.UnstructuredGrid(mri_mesh)

        # extract surface and restrict Eulerian mesh to whatever lies within it 
        mri_mesh_surface = mri_mesh_unstructured.extract_surface()


        selected = eulerian_mesh_centers.select_enclosed_points(mri_mesh_surface, tolerance=0.0, inside_out=False, check_surface=True)

        # all_scalars = True means all points must be within for a cell to stay 
        # all_scalars = False means any point must be within for a cell to stay 
        eulerian_mesh_centers_restricted = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=False) 


        # copy MRI data structure 
        eulerian_resampled = mri_mesh_unstructured.copy(deep=True)
        # add fields 
        eulerian_resampled.point_data['U'] = eulerian_resampled.point_data['velocity'] * 0.0
        eulerian_resampled.point_data['n_to_average'] = eulerian_resampled.point_data['image'] * 0    


        # use KD tree following example 
        # https://github.com/pyvista/pyvista-support/issues/107

        t0 = time.time()
        frac_to_process = 1

        kdtree_mri_mesh = KDTree(mri_mesh_unstructured.points.astype(np.double))

        # dist, mri_idx = kdtree_mri_mesh.query(eulerian_mesh_centers_restricted.points[0:math.floor(eulerian_mesh_centers.n_points/frac_to_process)])
        dist, mri_idx = kdtree_mri_mesh.query(eulerian_mesh_centers_restricted.points)

        print("mri_idx[0:40] ", mri_idx[0:40])

        # for sim_idx,pt in enumerate(eulerian_mesh_centers.points[0:math.floor(eulerian_mesh_centers.n_points/frac_to_process)] ):
        # for sim_idx in range(math.floor(eulerian_mesh_centers.n_points/frac_to_process)):
        for sim_idx in range(eulerian_mesh_centers_restricted.n_points):

            # mri_idx_temp = mri_idx[sim_idx]
            eulerian_resampled.point_data['U'][mri_idx[sim_idx]] += eulerian_mesh_centers_restricted.point_data['U'][sim_idx]
            eulerian_resampled.point_data['n_to_average'][mri_idx[sim_idx]] += 1 

        # compute the means 
        for i in range(eulerian_resampled.n_points):
            if eulerian_resampled.point_data['n_to_average'][i] != 0:
                eulerian_resampled.point_data['U'][i] /= eulerian_resampled.point_data['n_to_average'][i]


        t1 = time.time()
        total_time = t1-t0
        print("total_time = ", total_time)
        est_total_time_full = total_time * frac_to_process
        print("est_total_time_full = ", est_total_time_full)        

        print(eulerian_resampled.point_data['U'][0:20])
        print(eulerian_resampled.point_data['n_to_average'][0:20])

        eulerian_resampled.save('eulerian_mean_test.vtu')


        eulerian_resampled_read = pyvista.read('eulerian_mean_test.vtu')

        boundary_mesh_name = '2_1_HealthyGeometry_PlanarInletOutlet_smoothed_capped.stl'
        boundary_mesh = pyvista.read(boundary_mesh_name)
        
        # casting not allowed, need a surfbace mesh 
        # boundary_mesh = pyvista.UnstructuredGrid(boundary_mesh)

        selected = eulerian_resampled_read.select_enclosed_points(boundary_mesh, tolerance=0.0, inside_out=False, check_surface=True)

        # all_scalars = True means all points must be within for a cell to stay 
        # all_scalars = False means any point must be within for a cell to stay 
        eulerian_resampled_restricted = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

        eulerian_resampled_restricted.save('eulerian_mean_test_restricted.vtu')

        if np.isnan(eulerian_resampled_restricted.point_data['U']).any():
            print('NaN found in eulerian_resampled_restricted')






        # insanely slow! 
        # 
        # t0 = time.time()
        # frac_to_process = 1.0e5

        # for sim_idx,pt in enumerate(eulerian_mesh_centers.points[0:math.floor(eulerian_mesh_centers.n_points/frac_to_process)] ):

        #     mri_idx = mri_mesh_unstructured.find_closest_point(pt)

        #     eulerian_resampled.point_data['U'][mri_idx] += eulerian_mesh_centers.point_data['U'][sim_idx]

        #     eulerian_resampled.point_data['n_to_average'][mri_idx] += 1 

        # print(eulerian_resampled.point_data['U'][0:20])
        # print(eulerian_resampled.point_data['n_to_average'][0:20])

        # t1 = time.time()
        # total_time = t1-t0
        # print("total_time = ", total_time)
        # est_total_time_full = total_time * frac_to_process
        # print("est_total_time_full = ", est_total_time_full)    


if __name__ == '__main__':
    run_main()
