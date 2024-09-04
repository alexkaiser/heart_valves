from __future__ import print_function
import os, sys, subprocess, time


def fix_visit_files(viz_directory):
    ''' 
    Checks viz directory to find output increment and rewrites visit files 
    '''
                
    prev_dir = os.getcwd()
    os.chdir(viz_directory)
    
    # minimum nonzero step 
    min_nonzero = float('Inf')
    max_idx = 0
    
    for f in os.listdir('.'): 
        if f.startswith('visit_dump.'):
            step_num = int(f.split('.')[1])
            
            if (step_num is not 0):
                if step_num < min_nonzero: 
                    min_nonzero = step_num
                    
            if step_num > max_idx: 
                max_idx = step_num
    
    dumps = open('dumps.visit', 'w')
    lag   = open('lag_data.visit', 'w')

    state = 0
    stride = min_nonzero
    max_step = max_idx

    dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
    lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 

    for state in range(stride, max_step, stride):
        dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
        lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 

    state = max_step
    dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
    lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 

    dumps.close()
    lag.close()
    os.chdir(prev_dir)


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
                
                # clean up visit files to be consistent after restarts
                # if number_restarts > 0:
                fix_visit_files(viz_dir_name)

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

