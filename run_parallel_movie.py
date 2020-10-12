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



if __name__ == '__main__':

    print "call arguments: "
    for arg in sys.argv:
        print arg, " "


    if os.path.isfile('done.txt'):
        
        print 'done.txt found'
    
        for f in os.listdir('.'):
            if f.startswith('viz'): 
                
                print 'Found viz directory'
                
                os.chdir(f)
                
                viz_dir_name = os.getcwd()
                
                # clean up visit files to be consistent after restarts
                # if number_restarts > 0:
                fix_visit_files(viz_dir_name)

                if len(sys.argv) < 2:
                    raise InputError("Must provide a script name")
                session_file_name = sys.argv[1]

                if len(sys.argv) < 3:
                    raise InputError("Must specify n_procs")
                n_procs = int(sys.argv[2])

                if len(sys.argv) >= 4:
                    try:
                        view_clipping = float(sys.argv[3])
                    except ValueError:
                        print "Not a float"
                        view_clipping = None 
                else:
                    view_clipping = None

                call_str_base = 'visit -cli -nowin -s ~/copies_scripts/make_movie_generic.py ~/copies_scripts/'
                call_str_base += session_file_name + " "

                processes = []

                for i in range(n_procs):
                    call_str = call_str_base + ' ' + str(n_procs) + ' ' + str(i)

                    if view_clipping is not None: 
                        call_str += " " + str(view_clipping)

                    print "call str with n_procs = ", n_procs, " proc_num = ", i, "call_str = ", call_str

                    p = subprocess.Popen(call_str, shell=True)
                    print "started process, p = ", p

                    processes.append(p)

                print "out of call loop in orig python script "
                time.sleep(30)

                while True: 

                    print "new check"

                    all_done = True 
                    for p in processes:
                        if p.poll() is None:
                            print "process p = ", p, "returned None, still running"
                            all_done = False
                        else:
                            print "process p = ", p, "returned not None, is complete"

                    if all_done:
                        break

                    time.sleep(30)

    # get some output names
    cwd = os.getcwd()
    cwd_split = cwd.split('/')

    # name files after the job if easy
    # path is always /home/adk354/scratch/JOB_NAME
    if (len(cwd_split) >= 5):
        base_name = cwd_split[4]
    else:
        base_name = 'frames'

    # ffmpeg from here so variables are all in place
    movie_string = 'ffmpeg -framerate 60 -i '
    movie_string += base_name
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name + '.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print 'something wrong in movie make, call returned prematurely'


    # # reduce by 10x
    # # 60 input, 60 output is 10x slow motion
    # # 600 input, 60 output is real time 

    movie_string = 'ffmpeg -framerate 600 -i '
    movie_string += base_name
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name + '_real_time.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print 'something wrong in movie make, call returned prematurely'



