import sys
import os 
import shutil 
import pyvista
import numpy as np 
import multiprocessing

def add_faces(mesh_with_faces, file_name):

    base_name = file_name.rsplit('.', 1)[0]
    file_name_copy = base_name + "_orig_copy.vtu"
    shutil.copyfile(file_name, file_name_copy)

    mesh = pyvista.read(file_name)
    mesh_with_faces.points = mesh.points
    mesh_with_faces.save(file_name)


if __name__ == '__main__':

    process_pa = False 
    if process_pa:

        mesh_with_faces_name = '13_3_layer_dspt025_five_layer_inlet.vtu'

        if "_192_" in os.getcwd():
            mesh_with_faces_name = '12_three_layer_192.vtu'
        elif "_768_" in os.getcwd():
            mesh_with_faces_name = '12_three_layer_fix_extender_mesh_pt02_ext_pt014.vtu'

        base_name = 'vessel'

    # aorta default 
    else: 

        
        base_name = 'aorta_384'

        if os.path.isfile(base_name + '.vtu'):
            print("file found")
            mesh_with_faces_name = base_name + '.vtu'
        elif os.path.isfile('../' + base_name + '.vtu'):
            print("file found one dir")
            shutil.copy('../' + base_name + '.vtu', '.') 
            mesh_with_faces_name = base_name + '.vtu'
        else: 
            mesh_with_faces_name = "4_aorta_remeshed_pt25mm_3cm_extender_layers_constriction.vtu"

        if "_192_" in os.getcwd():
            mesh_with_faces_name = '6_aorta_remeshed_pt5mm_2cm_extender_layers_constriction.vtu'
            base_name = 'aorta_192'


    if not os.path.isfile(mesh_with_faces_name):
        if os.path.isfile(os.path.expanduser('~') + '/mitral_fully_discrete/' + mesh_with_faces_name):
            shutil.copy(os.path.expanduser('~') + '/mitral_fully_discrete/' + mesh_with_faces_name, '.') 
        else: 
            raise FileNotFoundError("cannot find mesh_with_faces_name file = ", mesh_with_faces_name)

    mesh_with_faces = pyvista.read(mesh_with_faces_name)



    use_multiprocessing = True
    if use_multiprocessing: 
        jobs = []

        for f in os.listdir('.'):
            if f.startswith(base_name) and f.endswith('.vtu'):
                if (not "_orig_copy" in f) and (mesh_with_faces_name not in f):
                    print("procssing file ", f)
                    
                    p = multiprocessing.Process(target=add_faces, args=(mesh_with_faces, f))
                    jobs.append(p)
                    p.start()
                    # add_faces(mesh_with_faces, f)

        for p in jobs:
            p.join()
    
    else:
        for f in os.listdir('.'):
            if f.startswith(base_name) and f.endswith('.vtu'):
                if not "_orig_copy" in f:
                    print("procssing file ", f)
                    add_faces(mesh_with_faces, f)

    # pool = multiprocessing.Pool() #use all available cores, otherwise specify the number you want as an argument
    # for i in range(nframes-1):add
    #     pool.apply_async(add_faces, args=(i,))
    # pool.close()
    # pool.join()
