import pyvista 
import os, sys, glob
import multiprocessing
import pdb

def write_pvd(base_name, dt, nsteps, extension, nprocs_sim=1):

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
                                     output_extension = ".vts",
                                     convert_to_cell_data=False,
                                     convert_to_point_data=False):

    print("processing frame_number, ", frame_number)

    eulerian_name = eulerian_basename_orig + str(frame_number).zfill(digits_output_eulerian) + eulerian_extension

    eulerian_mesh_orig = pyvista.read(eulerian_name)

    mri_name = mri_basename + str(frame_number).zfill(digits_output_mri) + mri_extension

    mri_mesh = pyvista.read(mri_name)

    eulerian_resampled = mri_mesh.sample(eulerian_mesh_orig)

    if convert_to_cell_data and convert_to_point_data:
        raise ValueError("Cannot convert_to_cell_data and convert_to_point_data")

    if convert_to_cell_data:
        eulerian_resampled = eulerian_resampled.point_data_to_cell_data()

    if convert_to_point_data:
        eulerian_resampled = eulerian_resampled.cell_data_to_point_data()
        suffix_eulerian += '_points'
    
    eulerian_name_out = eulerian_basename_orig + suffix_eulerian
    eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + output_extension

    eulerian_resampled.save(eulerian_name_out)




if __name__ == '__main__':

    run_all = False 

    if run_all:

        nframes = 20

        # jobs = []
        # for i in range(nframes):
        #     p = multiprocessing.Process(target=interpolate_eulerian_to_mri_mesh, args=(i,))
        #     jobs.append(p)
        #     p.start()

        # for p in jobs:
        #     p.join()

        eulerian_basename_orig="eulerian_vars_averaged"
        digits_output_eulerian = 4
        eulerian_extension = ".vtu"
        suffix_eulerian = '_resampled'
        mri_basename = "HealthyNative_ForSim_"
        digits_output_mri = 2
        mri_extension = ".vtk"
        output_extension = ".vts"
        convert_to_cell_data = False 
        convert_to_point_data = True

        pool = multiprocessing.Pool() #use all available cores, otherwise specify the number you want as an argument
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
                                                                     convert_to_cell_data,
                                                                     convert_to_point_data))
        pool.close()
        pool.join()


        cycle_duration = 8.3250000000000002e-01
        mri_read_times_per_cycle = 10 
        dt_mri_read = cycle_duration / mri_read_times_per_cycle
        output_times_per_cycle   = 20 
        dt_output = cycle_duration / output_times_per_cycle

        eulerian_basename_orig="eulerian_vars_averaged"

        suffix_eulerian = '_resampled'

        if convert_to_point_data: 
            suffix_eulerian += "_points"

        extension = "vts"

        write_pvd(eulerian_basename_orig + suffix_eulerian, dt_output, output_times_per_cycle, extension)



    tests = True
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



