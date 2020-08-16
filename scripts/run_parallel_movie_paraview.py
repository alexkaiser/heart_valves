
import os, sys, subprocess, time


def run_command_parallel(call_str_base, n_procs):

    processes = []

    for i in range(n_procs):
        call_str = call_str_base + ' ' + str(n_procs) + ' ' + str(i)

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



if __name__ == '__main__':

    print "call arguments: "
    for arg in sys.argv:
        print arg, " "

    script_dir = "~/copies_scripts/"
    # script_dir = "~/mitral_fully_discrete/scripts/"

    # get some output names
    cwd = os.getcwd()
    cwd_split = cwd.split('/')

    # name files after the job if easy
    # path is always /home/adk354/scratch/JOB_NAME
    if (len(cwd_split) >= 5):
        base_name = cwd_split[4]
    else:
        base_name = 'frames'
    base_name += "_paraview"

    if os.path.isfile('done.txt'):
        
        print 'done.txt found'
    
        for f in os.listdir('.'):
            if f.startswith('viz'): 
                
                print 'Found viz directory'
                
                os.chdir(f)
                
                viz_dir_name = os.getcwd()
                
                if len(sys.argv) < 2:
                    raise InputError("Must specify session_file_name")
                session_file_name = sys.argv[1]

                if len(sys.argv) < 3:
                    raise InputError("Must specify n_procs")
                n_procs = int(sys.argv[2])

                call_str_base = "pvpython " + session_file_name + " " + base_name + " "
                print "call_str_base = ", call_str_base 
                run_command_parallel(call_str_base, n_procs)


    # ffmpeg from here so variables are all in place
    movie_string = 'ffmpeg -y -framerate 60 -i '
    movie_string += base_name
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name + '.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print 'something wrong in movie make, call returned prematurely'


    # # reduce by 10x
    # # 60 input, 60 output is 10x slow motion
    # # 600 input, 60 output is real time 

    movie_string = 'ffmpeg -y -framerate 600 -i '
    movie_string += base_name
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name + '_real_time.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print 'something wrong in movie make, call returned prematurely'

