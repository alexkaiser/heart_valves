import pyvista 
import meshio
import os, sys, glob
import subprocess
import math 
from natsort import natsorted
import multiprocessing 
import pdb
import numpy as np 


def write_pvd(basename, dt, nsteps, extension, nprocs_sim=1):

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

                filename_out = basename + '.pvd'
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
                tmp_str += basename + str(n).zfill(4) + '/' # sorted into directories 
            tmp_str += basename + str(n).zfill(4) + '.' 
            if nprocs_sim > 1:
                tmp_str += str(proc) + '.'        
            tmp_str += extension
            tmp_str += '"/>\n'

            f_write.write(tmp_str)

    f_write.write(suffix)
    f_write.close()


def read_distributed_vtr(dir_name):

    files = natsorted(glob.glob(dir_name + "/*.vtr"))
    blocks = pyvista.MultiBlock([pyvista.RectilinearGrid(f) for f in files])
    return blocks.combine()


def sort_points_and_point_data(points, point_data=None):
    '''
    Sorts points by x,y,z
    Applies permutation to point and point_data
    '''

    indices = np.lexsort((points[:,0], points[:,1], points[:,2])) 

    points_sorted = points[indices]

    if point_data is not None:
        point_data_sorted = point_data[indices]
    else:
        point_data_sorted = None

    return points_sorted, point_data_sorted


def generate_cells(NX,NY,NZ):
    '''
    generates cells in meshio format 
    '''

    get_1d_idx = lambda i,j,k : i + j*NX + k*NX*NY

    cells = []

    # every point gets a cell 
    # no cells on the top row of points 
    for i in range(NX-1):
        for j in range(NY-1):
            for k in range(NZ-1):

                # vtk order points the normal of the "bottom" face up 
                # and the top face also up 
                cell_tmp = []
                cell_tmp.append(get_1d_idx(i  , j  , k  ))
                cell_tmp.append(get_1d_idx(i+1, j  , k  ))
                cell_tmp.append(get_1d_idx(i+1, j+1, k  ))
                cell_tmp.append(get_1d_idx(i  , j+1, k  ))
                cell_tmp.append(get_1d_idx(i  , j  , k+1))
                cell_tmp.append(get_1d_idx(i+1, j  , k+1))
                cell_tmp.append(get_1d_idx(i+1, j+1, k+1))
                cell_tmp.append(get_1d_idx(i  , j+1, k+1))

                cells.append(cell_tmp)

    return cells 



def convert_mesh_to_center_points(mesh, NX, NY, NZ):
    ''' 
    Takes cell data mesh and converts cell_array 'U' to point_array at center 
    '''

    mesh_cell_centers = mesh.cell_centers()
    points_sorted, point_data_sorted_U = sort_points_and_point_data(mesh_cell_centers.points, mesh_cell_centers.point_arrays['U'])
    cells_sorted = generate_cells(NX,NY,NZ)

    mesh_points_meshio_format = meshio.Mesh(points_sorted, 
                                        cells={"hexahedron": cells_sorted}, 
                                        point_data={'U': point_data_sorted_U})

    return pyvista.wrap(mesh_points_meshio_format)



def combine_mesh(basename,                           
                 range_to_output, 
                 label='_cells', 
                 extension='vtu', 
                 point_data=False, 
                 NX=None, 
                 NY=None, 
                 NZ=None, 
                 proc_num=0, 
                 nprocs=1):

    if isinstance(range_to_output, int):
        to_process = range(range_to_output)
    else:
        to_process = range_to_output

    if point_data:
        if (NX is None) or (NY is None) or (NZ is None):
            raise ValueError("Must provide values for NX,NY,NZ when calling points")

    for i in to_process:
        if (i % nprocs) == proc_num:

            dir_name = basename + str(i).zfill(4)

            fname_out = basename + label + str(i).zfill(4) + '.' + extension

            # read distributed vtr 
            mesh = read_distributed_vtr(dir_name)

            if point_data:
                mesh = convert_mesh_to_center_points(mesh, NX, NY, NZ)
                
            mesh.save(fname_out)


def read_and_convert_points(basename,
                            frame_number, 
                            NX, 
                            NY,
                            NZ,
                            label='_points', 
                            extension='vtu'):

    if point_data:
        if (NX is None) or (NY is None) or (NZ is None):
            raise ValueError("Must provide values for NX,NY,NZ when calling points")

    fname_in = basename + str(frame_number).zfill(4) + '.' + extension

    fname_out = basename + label + str(frame_number).zfill(4) + '.' + extension

    mesh = pyvista.read(fname_in)

    if point_data:
        mesh = convert_mesh_to_center_points(mesh, NX, NY, NZ)
        
    mesh.save(fname_out)




if __name__ == '__main__':

    if len(sys.argv) >= 2:
        nprocs = int(sys.argv[1]) # number of threads to launch  
    else: 
        print("using default nprocs = 1")
        nprocs = 1 

    # crop times for debug 
    # times = times[:5]

    point_data = True 


    coarse = True
    med_coarse = False
    fine = False
    xfine = False
    xxfine = False

    if coarse: 
        NX = 104
        NY = 72
        NZ = 64
        frame_number_temp = 1393
    elif med_coarse:
        NX = 138
        NY = 96
        NZ = 85
        frame_number_temp = 1370
    elif fine: 
        NX=208
        NY=144
        NZ=128
        frame_number_temp = 1350
    elif xfine: 
        NX = 277
        NY = 192
        NZ = 170
        frame_number_temp = 1350 # or 1351 
    elif xxfine:
        NX=332
        NY=230
        NZ=204
        frame_number_temp = 1350 # or 1351 
    else:
        raise ValueError('must provide values')

    basename = "eulerian_vars"
    
    if point_data:
        label = '_points'
    else: 
        label = '_cells'

    extension = 'vtu'

    run_all = False 
    if run_all:

        # first make sure there is a times file 
        if not os.path.isfile('times.txt'):
            subprocess.call('visit -cli -nowin -s ~/heart_valves/scripts/write_times_file_visit.py', shell=True)

        times = []
        times_file = open('times.txt', 'r')
        for line in times_file:
            times.append(float(line)) 

        if times[0] == 0:
            dt = times[1]
        else: 
            dt = times[1] - times[0]

        nsteps = len(times)

        jobs = []
        for proc_num in range(nprocs):
            p = multiprocessing.Process(target=combine_mesh, args=(basename, 
                                                                   nsteps, 
                                                                   label, 
                                                                   extension,
                                                                   point_data, 
                                                                   NX,
                                                                   NY, 
                                                                   NZ, 
                                                                   proc_num, 
                                                                   nprocs))
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

        basename_out = basename + label

        write_pvd(basename_out, dt, nsteps, extension, nprocs_sim=1)


    run_convert_only = False
    if run_convert_only:

        cycle_duration = 8.3250000000000002e-01
        mri_read_times_per_cycle = 10 
        dt_mri_read = cycle_duration / mri_read_times_per_cycle

        nsteps = 10

        basename += '_averaged'

        use_multi = True 
        if use_multi: 
            jobs = []
            for frame_number in range(nsteps):

                p = multiprocessing.Process(target=read_and_convert_points, args=(basename,
                                                                                  frame_number, 
                                                                                  NX,
                                                                                  NY,
                                                                                  NZ,
                                                                                  label,
                                                                                  extension))
                jobs.append(p)
                p.start()

            for p in jobs:
                p.join()
        else: 
            for frame_number in range(nsteps):
            # for frame_number in [0]:
                read_and_convert_points(basename,
                                        frame_number, 
                                        NX, 
                                        NY,
                                        NZ,
                                        label, 
                                        extension)

        basename_out = basename + label

        write_pvd(basename_out, dt_mri_read, nsteps, extension, nprocs_sim=1)



    local = True 
    if local:

        # range_to_output = [0]
        # range_to_output = [1350]
        range_to_output = [0, frame_number_temp]
        combine_mesh(basename, 
                     range_to_output, 
                     label, 
                     extension,
                     point_data, 
                     NX,
                     NY, 
                     NZ)



