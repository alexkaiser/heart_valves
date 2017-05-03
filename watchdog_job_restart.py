


import subprocess  
import time 
import sys 
import os

def full_run_line(run_line, input_name, options, restart_number=None, restart_dir=None): 
    ''' 
    Returns the string to run the given code 
    
    If restart_number = 0, then no restart string is included 
    
    If restart_dir is None, then the directory is called restart_IB3d_tree_cycle
    '''
    
    line = run_line + ' ' + input_name
    
    if (restart_number > 0) and (restart_number is not None):
        if restart_dir is None: 
            dir_name = 'restart_IB3d_tree_cycle'
        else: 
            dir_name = restart_dir
        line += ' ' + dir_name + ' ' + str(restart_number)
        
    line += ' ' + options
    return line


def log_file_exists(process=None):
    ''' check in loop to make sure log file is found, sleep for a minute 100 times if not'''
    
    for i in range(10):
        if os.path.isfile('IB3d.log'):
            print 'log file found, script moving forward'
            return True
        
        if process is not None:
            if process.poll() is not None:
                print 'process reported stopped, returning from log_file_exists'
                return False
        
        
        time.sleep(60)

    return False


def get_restart_number(restart_dir):
    ''' Gets the number of the most recent restart from the restart_dir '''

    if not os.path.exists(restart_dir):
        return None

    largest_restart_num = None
    for f in os.listdir(restart_dir):
        str_list = f.split('.')
        if str_list[0] == 'restore':
            val = int(str_list[1])
            if largest_restart_num is None:
                largest_restart_num = val
            if val > largest_restart_num:
                largest_restart_num = val

    if largest_restart_num is None:
        print 'No restart files found, restart at beginning'

    return largest_restart_num




if __name__ == '__main__':
    
    assert (len(sys.argv) >= 4) 
        
    run_line    = sys.argv[1]
    input_name  = sys.argv[2]
    options     = sys.argv[3]    
    
    print 'run line = ', run_line 
    print 'input_name =', input_name
    print 'options = ', options  
    
    if len(sys.argv) >= 5:
        restart_dir = sys.argv[4]
    else:
        print ''
        restart_dir = 'restart_IB3d_tree_cycle'
        
    # clean up the old done file if needed 
    if os.path.isfile('done.txt'):
        code = subprocess.call('rm done.txt', shell=True
            if code is None:
                print 'removal of done.txt failed\n'

    # check if we have restart available,
    # returns None if not
    restart_number = get_restart_number(restart_dir)
    
    to_run = full_run_line(run_line, input_name, options, restart_number, restart_dir)
    
    print 'full run line = ', to_run
    
    script_name = 'run_script_first.sh'
    script = open(script_name, 'w')
    
    
    script_prelims = '''
    #!/bin/bash
    
    env_log=env.log
    while read -r line; do
        export "$line"
    done < $env_log
        
    '''
    
#    script_prelims_restart = '''
#    unset PBS_JOBID
#    export PBS_JOBID=
    
#    '''

    script.write(script_prelims)
    script.write(to_run + '\n\n')
    script.close()
    
    print 'python thinks current working dir is ', os.getcwd()
    
    script_call_line = 'sh ' + script_name
    print 'Popen will be called with ', script_call_line
    
    # start the orig process
    print 'to process...'
    current_sh_calls_mpi = subprocess.Popen(script_call_line, shell=True)
    print 'should be running, control back to python'

    # check for file
    log_found = log_file_exists()
    if not log_found:
        print 'No log file found after 10 checks, killing python script'
        sys.exit()

    wait_time_s = 10*60              # check every ten minutes, max 20 minutes lost
    wait_time_before_restart = 2*60  # after everything is killed, just hang out for two minutes
    wait_time_after_restart = 5*60   # once the restart goes, add a few extra minutes for initialization 
    number_restarts = 0
    prev_time = os.path.getmtime('IB3d.log')
    check_number = 0

    # wait, check if stopped, restart if needed
    while True:
        
        print 'On loop check number ', check_number
        
        # just hang out
        time.sleep(wait_time_s)
        
        # last changes to log
        mod_time = os.path.getmtime('IB3d.log')
        
        # check whether the job has completed
        # no restart should occur if so
        # poll returns none if still running
        if current_sh_calls_mpi.poll() is not None:
            print 'MPI run has stopped, check for completion or crashes.'
            print 'Killing python watchdog.'
            sys.exit()
        
        if mod_time == prev_time:
            print 'On check number ', check_number, ', modification time unchanged'
            print 'Initiating restart.'

            # kill the sh script
            current_sh_calls_mpi.kill()
            current_sh_calls_mpi.wait()
            if current_sh_calls_mpi.poll() is None:
                print 'Shell script ', script_name, ' is still running (even though it should have completed).'
                print 'Beware of strange behavior'

            # kill the MPI jobs with a shell script, wait for this to finish
            code = subprocess.call('sh kill_all_mpi.sh', shell=True)
            if code is None:
                print 'kill_all_mpi is still running (even though it should have completed).'
                print 'Beware of strange behavior'

            # move old log and output files
            move_str = 'mv IB3d.log IB3d.log_restart_' + str(number_restarts)
            subprocess.call(move_str, shell=True)
            move_str = 'mv output.txt output.txt_restart_' + str(number_restarts)
            subprocess.call(move_str, shell=True)

            # make a new script
            script_name = 'run_script_restart_' + str(number_restarts) + '.sh'
            script = open(script_name, 'w')
            script.write(script_prelims)
            
            restart_number = get_restart_number(restart_dir)
            print 'restart_number = ', restart_number

            to_run = full_run_line(run_line, input_name, options, restart_number, restart_dir)
            script.write(to_run + '\n\n')
            script.close()

            time.sleep(wait_time_before_restart)
            
            # kill the MPI jobs with a shell script (again...), wait for this to finish
            code = subprocess.call('sh kill_all_mpi.sh', shell=True)
            if code is None:
                print 'kill_all_mpi is still running (even though it should have completed).'
                print 'Beware of strange behavior'
            time.sleep(wait_time_before_restart)
            
            
            print 'restart number ', number_restarts, 'at step ', number_restarts

            
            current_sh_calls_mpi = subprocess.Popen('sh ' + script_name, shell=True)

            # make sure we have the log file again
            log_found = log_file_exists()
            if not log_found:
                print 'No log file found after 10 checks, try to submit again'
                    
                current_sh_calls_mpi.kill()
                current_sh_calls_mpi.wait()
                if current_sh_calls_mpi.poll() is None:
                    print 'Shell script ', script_name, ' is still running (even though it should have completed).'
                    print 'Beware of strange behavior'
                
                time.sleep(wait_time_before_restart)
                
                # try again
                current_sh_calls_mpi = subprocess.Popen('sh ' + script_name, shell=True)
                
                log_found_again = log_file_exists()
                if not log_found_again:
                    print 'No log file found after 10 checks, cancel job.'
                    sys.exit()

            # extra wait to allow for initialization
            time.sleep(wait_time_after_restart)

            number_restarts += 1


        prev_time = mod_time
        check_number += 1
    
    # submit movie script for post processing 
    if os.path.isfile('done.txt'):
        for f in os.listdir('.'):
            if f.startswith('viz'): 
                
                os.chdir(f)
                
                viz_dir_name = os.getcwd()
                
                movie_script = open('make_movie.sbatch', 'w')
                
                slurm = '''#!/bin/bash
                #SBATCH --nodes=1
                #SBATCH --ntasks=1
                #SBATCH --time=4:00:00
                #SBATCH --mem=16GB
                #SBATCH --job-name=movie_post_process
                #SBATCH --mail-user=kaiser@cims.nyu.edu
                #SBATCH --mail-type=ALL
                
                ''' 
                
                movie_script.write(slurm)
                
                movie_script.write('\n')                
                movie_script.write('cd ' + viz_dir_name + ' \n')
                movie_script.write('visit -cli -nowin -s ~/mitral_fully_discrete/make_three_slice_movie.py \n')          
                
                movie_script.close()
                
                code = subprocess.call('sbatch make_movie.sbatch', shell=True)
                if code is None:
                    print 'submit of movie script failed, check for problems.\n'
                
                break 
    
    print 'done with main'


