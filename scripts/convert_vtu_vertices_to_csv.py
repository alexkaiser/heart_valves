import pyvista 
import os
import numpy as np 
import subprocess
import multiprocessing


def convert_csv(basename, frame_number, extension_in='.vtu', extension_out='.csv'):

    fname_in = basename + str(frame_number).zfill(4) + extension_in
    fname_out = basename + str(frame_number).zfill(4) + extension_out
    
    mesh = pyvista.read(fname_in)

    points = mesh.points 

    np.savetxt(fname_out, mesh.points, delimiter=', ')



if __name__ == '__main__':


    # basename = "aortic_no_partition_"

    # if '_192_' in os.getcwd(): 
    #     resolution_string = '192'
    # if '_384_' in os.getcwd(): 
    #     resolution_string = '384'

    # basename += resolution_string

    # first make sure there is a times file 
    if not os.path.isfile('times.txt'):
        subprocess.call('visit -cli -nowin -s ~/heart_valves/scripts/write_times_file_visit.py', shell=True)

    times = []
    times_file = open('times.txt', 'r')
    for line in times_file:
        times.append(float(line)) 
    nsteps = len(times)

    frame_number = 444 

    lag_name_base_to_check = ['aortic']

    for lag_file in os.listdir('..'):
        for lag_base in lag_name_base_to_check: 
            if lag_file.startswith(lag_base) and lag_file.endswith('.vertex'):

                basename = lag_file.rsplit('.', 1)[0]
                # convert_csv(basename, frame_number)

                pool = multiprocessing.Pool() #use all available cores, otherwise specify the number you want as an argument
                for i in range(nsteps):
                    pool.apply_async(convert_csv, args=(basename, i))
                pool.close()
                pool.join()








