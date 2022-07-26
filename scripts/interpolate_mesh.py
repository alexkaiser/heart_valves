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

    run_all = True 

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

        mri_name = mri_basename + str(frame_number).zfill(digits_output_mri) + mri_extension

        mri_mesh = pyvista.read(mri_name)

        eulerian_resampled = mri_mesh.sample(eulerian_mesh_orig)

        eulerian_resampled = eulerian_resampled.point_data_to_cell_data()

        suffix_eulerian = '_resampled'
        eulerian_name_out_extension = '.vtk'
        eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension

        eulerian_resampled.save(eulerian_name_out)

        # gaussian at scan radius 
        # does not work written like this, unclear why 
        # radius = 0.09 * 4
        # eulerian_mesh_guass_interp = mri_mesh.interpolate(eulerian_mesh_orig, radius=radius)
        # suffix_eulerian = '_gaussian_interp'
        # eulerian_name_out = eulerian_basename_orig + suffix_eulerian
        # eulerian_name_out += str(frame_number).zfill(digits_output_eulerian) + eulerian_name_out_extension
        # eulerian_mesh_guass_interp.save(eulerian_name_out)



