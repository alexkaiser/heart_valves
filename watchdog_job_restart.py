


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
    
    str = run_line + ' ' + input_name
    
    if restart_number > 0: 
        if restart_dir is None: 
            dir_name = 'restart_IB3d_tree_cycle'
        else: 
            dir_name = restart_dir
        str += ' ' + dir_name + ' ' + str(restart_number)
        
    str += ' ' + options
    return str
    
    


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
    
    to_run = full_run_line(run_line, input_name, options)
    
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
    
    script.write(script_prelims)
    script.write(to_run + '\n\n') 
    
    print 'python thinks current working dir is ', os.getcwd()
    
    script_call_line = 'sh ' + script_name
    print 'Popen will be called with ', script_call_line
    
    # start the orig process
    print 'to process...'
    subprocess.Popen(script_call_line, shell=True)
    print 'should be running, control back to python'
    
    # wait, check if stopped, restart if needed
    while true:
    
    
    
    
    
    
    
    print 'done with main'






    



