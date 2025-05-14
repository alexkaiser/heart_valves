import sys, os
import subprocess


def write_pvd(basename, times, nsteps, extension, nprocs_sim=1, skip_last_step=True, file_suffix=""):

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

                filename_out = basename + file_suffix + '.pvd'
                print("filename_out = ", filename_out)

                f_write = open(filename_out, 'w')
                f_write.write(prefix)
                initialized = True

            tmp_str = '    <DataSet timestep="'
            tmp_str += '{:.14f}'.format(times[n])
            tmp_str += '" group="" part="'
            tmp_str += str(proc) + '"'

            tmp_str += ' file="'
            if nprocs_sim > 1:
                tmp_str += basename + str(n).zfill(4) + '/' # sorted into directories 
            tmp_str += basename + str(n).zfill(4) + file_suffix + '.' 
            if nprocs_sim > 1:
                tmp_str += str(proc) + '.'        
            tmp_str += extension
            tmp_str += '"/>\n'

            f_write.write(tmp_str)

    f_write.write(suffix)
    f_write.close()



if __name__ == '__main__':

    # first make sure there is a times file 
    if not os.path.isfile('times.txt'):
        subprocess.call('visit -cli -nowin -s ~/copies_scripts/write_times_file_visit.py', shell=True)

    times = []
    times_file = open('times.txt', 'r')
    for line in times_file:
        times.append(float(line)) 
    
    nsteps = len(times)

    # names = ['eulerian_vars_restricted_points', 'aorta_384', 'aortic_no_partition_384', 'aortic_no_partition_384_cylinder']
    # names = ['aortic_no_partition_384', 'aortic_no_partition_384_cylinder']

    lag_name_base_to_check = ['aortic', 'aorta']

    for lag_file in os.listdir('..'):
        if lag_file.endswith('.vertex'):
            for lag_base in lag_name_base_to_check: 
                if lag_file.startswith(lag_base):

                    basename = lag_file.rsplit('.', 1)[0]
                    file_suffix_list = ["", "_faces"]
                    extension = 'vtu'

                    for file_suffix in file_suffix_list: 
                        print("writing basename = ", basename, "file_suffix = ", file_suffix)
                        write_pvd(basename, times, nsteps, extension, nprocs_sim=1, skip_last_step=False, file_suffix=file_suffix)

    if "fish" not in os.getcwd(): 
        if "graft" not in os.getcwd(): 
            basename = 'eulerian_vars_restricted_points'
            write_pvd(basename, times, nsteps, extension, nprocs_sim=1, skip_last_step=False)            



    basename = 'eulerian_vars'
    extension = 'vtr'

    if '_192_' in os.getcwd(): 
        nprocs_sim = 24
    elif '_384_' in os.getcwd(): 
        nprocs_sim = 48
    else: 
        print("Using default nprocs_sim of 48")
        nprocs_sim = 48

    print("writing eulerian_vars with nprocs_sim = ", nprocs_sim)

    write_pvd(basename, times, nsteps, extension, nprocs_sim=nprocs_sim, skip_last_step=False)

    # basename = 'aortic_no_partition_384'
    # file_suffix = '_faces'

    # write_pvd(basename, times, nsteps, extension, nprocs_sim=1, skip_last_step=False, file_suffix=file_suffix)



