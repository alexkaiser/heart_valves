from __future__ import print_function
import os, sys, subprocess, time


def run_command_parallel(call_str_base, n_procs):

    processes = []

    for i in range(n_procs):
        call_str = call_str_base + ' ' + str(n_procs) + ' ' + str(i)

        print("call str with n_procs = ", n_procs, " proc_num = ", i, "call_str = ", call_str)

        p = subprocess.Popen(call_str, shell=True)
        print("started process, p = ", p)

        processes.append(p)

    print("out of call loop in orig python script ")
    time.sleep(30)

    while True: 

        print("new check")

        all_done = True 
        for p in processes:
            if p.poll() is None:
                print("process p = ", p, "returned None, still running")
                all_done = False
            else:
                print("process p = ", p, "returned not None, is complete")

        if all_done:
            break

        time.sleep(30)


if __name__ == '__main__':

    print("call arguments: ")
    for arg in sys.argv:
        print(arg, " ")

    script_dir = "~/copies_scripts/"
    # script_dir = "~/mitral_fully_discrete/scripts/"

    lag_name_base_to_check = ['aortic', 'vessel', 'aorta_384.', 'aorta_192.']

    if os.path.isfile('done.txt'):
        
        print('done.txt found')
    
        for f in os.listdir('.'):
            if f.startswith('viz'): 
                
                print('Found viz directory')
                
                os.chdir(f)
                
                viz_dir_name = os.getcwd()

                if len(sys.argv) < 2:
                    raise InputError("Must specify n_procs_sim")
                n_procs_sim = int(sys.argv[1])

                if len(sys.argv) < 3:
                    raise InputError("Must specify n_procs")
                n_procs = int(sys.argv[2])

                call_str_base = 'visit -cli -nowin -s ' + script_dir + 'export_eulerian_visit_to_vtk.py '
                call_str_base += " eulerian_vars vtr " + str(n_procs_sim) + " "

                run_command_parallel(call_str_base, n_procs)

                for lag_file in os.listdir('..'):
                    if lag_file.endswith('.vertex'):
                        for lag_base in lag_name_base_to_check: 
                            if lag_file.startswith(lag_base):

                                print("found lag file ", lag_file, ", processing parallel")
                                base_name_lag = lag_file.rsplit('.', 1)[0]

                                print("base_name_lag = ", base_name_lag)

                                call_string_lag = 'visit -cli -nowin -s ' + script_dir + 'export_lag_visit_to_vtk.py '
                                call_string_lag += base_name_lag 
                                call_string_lag += " vtu " # always vtu for silo files 

                                run_command_parallel(call_string_lag, n_procs)                        

