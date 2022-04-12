import pyvista 
import os, sys, glob
import subprocess
import math 
from natsort import natsorted
import multiprocessing 

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



def read_distributed_vtr(dir_name):

    files = natsorted(glob.glob(dir_name + "/*.vtr"))
    #print("files = ", files)
    blocks = pyvista.MultiBlock([pyvista.read(f) for f in files])
    return blocks.combine()


def average_eulerian_mesh_one_step(idx_mri_read, eulerian_var_names, times, cycle_duration, cycles_to_output, dt_mri_read, base_dir, base_name_out, extension):

    # always start with initial eulerian mesh 
    dir_name = "eulerian_vars" + str(0).zfill(4)

    mesh = read_distributed_vtr(dir_name)
    n_to_average = 0
    for var_name in eulerian_var_names: 
        mesh[var_name] *= 0.0

    # average over times 
    for idx, t in enumerate(times): 

        # check if time in range 
        cycle_num = math.floor(t / cycle_duration)

        # skip cycle one 
        if cycle_num in cycles_to_output:
            
            dir_name = "eulerian_vars" + str(idx).zfill(4)

            # time since start of this cycle 
            t_reduced = t % cycle_duration

            idx_mri_read_temp = math.floor(t_reduced / dt_mri_read)

            if idx_mri_read == idx_mri_read_temp:
                print("processing step ", idx)

                mesh_tmp = read_distributed_vtr(dir_name)

                for var_name in eulerian_var_names: 
                    mesh[var_name] += mesh_tmp[var_name]
                    
                n_to_average += 1.0

            # print("t = ", t, "t_reduced = ", t_reduced, "idx_mri_read = ", idx_mri_read)

    for var_name in eulerian_var_names:
        # if none to average, mesh output as all zeros  
        if n_to_average != 0:
            mesh[var_name] /= float(n_to_average)

    fname = base_name_out + str(idx_mri_read).zfill(4) + '.' + extension
    mesh.save(base_dir + "/" + fname)    



if __name__ == '__main__':

    if len(sys.argv) >= 2:
        nprocs_sim = int(sys.argv[1]) # number of procs in the sim, which determines how many files go into the decomposed data
    else: 
        print("using default nprocs_sim = 1")
        nprocs_sim = 1 

    # first make sure there is a times file 
    if not os.path.isfile('times.txt'):
        subprocess.call('visit -cli -nowin -s ~/copies_scripts/write_times_file_visit.py', shell=True)

    times = []
    times_file = open('times.txt', 'r')
    for line in times_file:
        times.append(float(line)) 

    eulerian = True
    lagrangian = True 

    first_cycle = False
    second_cycle = True
    third_cycle = False
    second_third_cycle = False

    if first_cycle:
        cycles_to_output = [0] # zero indexed 
        # set up some directories 
        base_dir = "vis_data_averaged_cycle_1"
    elif second_cycle:
        cycles_to_output = [1] # zero indexed 
        # set up some directories 
        base_dir = "vis_data_averaged_cycle_2"
    elif third_cycle:
        cycles_to_output = [2] # zero indexed 
        # set up some directories 
        base_dir = "vis_data_averaged_cycle_3"
    elif second_third_cycle:
        cycles_to_output = [1,2] # zero indexed 
        # set up some directories 
        base_dir = "vis_data_averaged_cycle_2_3"
    else:
        cycles_to_output = [1,2,3] # zero indexed 
        # set up some directories 
        base_dir = "vis_data_averaged_cycle_2_3_4"


    cycle_duration = 8.3250000000000002e-01
    mri_read_times_per_cycle = 10 
    dt_mri_read = cycle_duration / mri_read_times_per_cycle
    output_times_per_cycle   = 20 
    dt_output = cycle_duration / output_times_per_cycle



    if not os.path.exists(base_dir):
        os.mkdir(base_dir)

    if eulerian:

        eulerian_var_names = ['P','Omega', 'U']

        # output file extension 
        extension = 'vtu'

        suffix = "_averaged"

        base_name_out = "eulerian_vars_mri_freq"

        # average all the Eulerian files here 
        # for idx_mri_read in range(mri_read_times_per_cycle):
        #     average_eulerian_mesh_one_step(idx_mri_read, eulerian_var_names, times, cycle_duration, cycles_to_output, dt_mri_read, base_dir, base_name_out, extension)

        jobs = []
        for idx_mri_read in range(mri_read_times_per_cycle):
            p = multiprocessing.Process(target=average_eulerian_mesh_one_step, args=(idx_mri_read, eulerian_var_names, times, cycle_duration, cycles_to_output, dt_mri_read, base_dir, base_name_out, extension))
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()


        # for idx_output in range(output_times_per_cycle):
        #     eulerian_dir_name = base_dir + '/' + 'eulerian_vars' + suffix + str(idx_output).zfill(4)
        #     if not os.path.exists(eulerian_dir_name):
        #         os.mkdir(eulerian_dir_name)

        # only average cycle 2 
        # cycles_to_include = [2] 

        # loops over parallel data structure as outer loop 
        # for proc_num in range(nprocs_sim):

        # read and zero meshes to use to accumulate from first mesh 
        dir_name = "eulerian_vars" + str(0).zfill(4)


        # read all time zero meshes 
        # meshes_mri_read = []
        # n_to_average = []
        # for idx_mri_read in range(mri_read_times_per_cycle):
        #     meshes_mri_read.append(read_distributed_vtr(dir_name))
        #     n_to_average.append(0)
        #     for var_name in eulerian_var_names: 
        #         meshes_mri_read[idx_mri_read][var_name] *= 0.0

        meshes_mri_read = []
        for idx_mri_read in range(mri_read_times_per_cycle):
            fname = base_name_out + str(idx_mri_read).zfill(4) + '.' + extension
            meshes_mri_read.append( pyvista.read(base_dir + "/" + fname) )

        meshes_output = []
        for idx_output in range(output_times_per_cycle):
            meshes_output.append(read_distributed_vtr(dir_name))
            for var_name in eulerian_var_names: 
                meshes_output[idx_output][var_name] *= 0.0

        # # average over times 
        # for idx, t in enumerate(times): 

        #     # check if time in range 
        #     cycle_num = math.floor(t / cycle_duration)

        #     # skip cycle one 
        #     if cycle_num in cycles_to_output:
        #         print("processing step ", idx)

        #         dir_name = "eulerian_vars" + str(idx).zfill(4)

        #         # time since start of this cycle 
        #         t_reduced = t % cycle_duration

        #         idx_mri_read = math.floor(t_reduced / dt_mri_read)

        #         mesh_tmp = read_distributed_vtr(dir_name)

        #         for var_name in eulerian_var_names: 
        #             meshes_mri_read[idx_mri_read][var_name] += mesh_tmp[var_name]
                    
        #         n_to_average[idx_mri_read] += 1.0

        #         # print("t = ", t, "t_reduced = ", t_reduced, "idx_mri_read = ", idx_mri_read)

        # print("n_to_average = ", n_to_average)

        # # convert sums to averages
        # for idx_mri_read in range(mri_read_times_per_cycle):
        #     for var_name in eulerian_var_names: 
        #         meshes_mri_read[idx_mri_read][var_name] /= float(n_to_average[idx_mri_read])

        # linearly interpolate before output         
        for idx_mri_read in range(mri_read_times_per_cycle):
            for var_name in eulerian_var_names: 
                meshes_output[2*idx_mri_read][var_name] = meshes_mri_read[idx_mri_read][var_name]

        for idx_mri_read in range(mri_read_times_per_cycle):
            idx_mri_read_next = (idx_mri_read + 1) % mri_read_times_per_cycle
            for var_name in eulerian_var_names: 
                meshes_output[2*idx_mri_read + 1][var_name] = 0.5 * (meshes_mri_read[idx_mri_read][var_name] + meshes_mri_read[idx_mri_read_next][var_name])

        for idx_output in range(output_times_per_cycle):
            eulerian_dir_name = base_dir 
            fname = "eulerian_vars" + suffix + str(idx_output).zfill(4) + '.' + extension
            meshes_output[idx_output].save(eulerian_dir_name + "/" + fname)

        # summary file         
        nprocs_output = 1
        write_pvd("eulerian_vars" + suffix, dt_output, output_times_per_cycle, extension, nprocs_output)
        os.rename("eulerian_vars" + suffix + '.pvd', base_dir + "/eulerian_vars" + suffix + '.pvd')


    if lagrangian:
        suffix = "_averaged"

        for lag_file in os.listdir('..'):
            if lag_file.endswith('.vertex'):

                print("found lag file ", lag_file, ", processing ")
                base_name_lag = lag_file.rsplit('.', 1)[0]
                print("base_name_lag = ", base_name_lag)

                # read and zero meshes to use to accumulate from first mesh 
                fname    = base_name_lag + str(0).zfill(4) + '.vtu'

                if not os.path.isfile(fname):
                    print("vtu file not found, cannot process this file, continuing")
                    continue 

                meshes_mri_read = []
                n_to_average = []
                for idx_mri_read in range(mri_read_times_per_cycle):
                    meshes_mri_read.append(pyvista.read(fname))
                    n_to_average.append(0)
                    meshes_mri_read[idx_mri_read].points *= 0.0

                meshes_output = []
                for idx_output in range(output_times_per_cycle):
                    meshes_output.append(pyvista.read(fname))
                    meshes_output[idx_output].points *= 0.0

                # average over times 
                for idx, t in enumerate(times): 

                    # check if time in range 
                    cycle_num = math.floor(t / cycle_duration)

                    # skip cycle one 
                    if cycle_num in cycles_to_output:

                        fname    = base_name_lag + str(idx).zfill(4) + '.vtu'

                        # time since start of this cycle 
                        t_reduced = t % cycle_duration

                        idx_mri_read = math.floor(t_reduced / dt_mri_read)

                        mesh_tmp = pyvista.read(fname)

                        meshes_mri_read[idx_mri_read].points += mesh_tmp.points
                            
                        n_to_average[idx_mri_read] += 1.0

                        # print("t = ", t, "t_reduced = ", t_reduced, "idx_mri_read = ", idx_mri_read)

                print("n_to_average = ", n_to_average)

                # convert sums to averages
                for idx_mri_read in range(mri_read_times_per_cycle):
                    meshes_mri_read[idx_mri_read].points /= float(n_to_average[idx_mri_read])

                # linearly interpolate before output
                for idx_mri_read in range(mri_read_times_per_cycle):
                    meshes_output[2*idx_mri_read].points = meshes_mri_read[idx_mri_read].points

                for idx_mri_read in range(mri_read_times_per_cycle):
                    idx_mri_read_next = (idx_mri_read + 1) % mri_read_times_per_cycle
                    meshes_output[2*idx_mri_read + 1].points = 0.5 * (meshes_mri_read[idx_mri_read].points + meshes_mri_read[idx_mri_read_next].points)

                for idx_output in range(output_times_per_cycle):
                    fname = base_name_lag + suffix + str(idx_output).zfill(4) + '.vtu'
                    meshes_output[idx_output].save(base_dir + "/" + fname)
                    # os.rename(fname, base_dir + "/" + base_name_lag + suffix + '.pvd')

                # summary file 
                extension = 'vtu'
                write_pvd(base_name_lag + suffix, dt_output, output_times_per_cycle, extension, 1)
                os.rename(base_name_lag + suffix + '.pvd', base_dir + "/" + base_name_lag + suffix + '.pvd')


























