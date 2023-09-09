import subprocess 
import sys 
import os 


if __name__ == '__main__':

    print "call arguments: "
    for arg in sys.argv:
        print(arg, " ")

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
        
        print('done.txt found')
    
        for f in os.listdir('.'):
            if f.startswith('viz'): 
                
                print('Found viz directory')
                
                os.chdir(f)
                
                viz_dir_name = os.getcwd()

                code = subprocess.call('cp ~/copies_scripts/run_aortic_plots.m .', shell=True)
                if code is None:
                    print('call returned prematurely')
                
                stride_found = False 
                if os.path.isfile('dumps.visit'):
                    visit_file = open('dumps.visit', 'r')
                    # first has zeros 
                    visit_file.readline()
                    visit_str = visit_file.readline()
                    visit_str = visit_str.split('.')[1] # format visit_dump.00000/summary.samrai 
                    stride = visit_str.split('/')[0]                    
                    print("stride = ", stride)
                    stride_found = True 

                call_str_matlab = 'matlab -nodisplay -nodesktop -r "run_aortic_plots('
                if 'high_pressure' in base_name:
                    call_str_matlab += "'high'" # outer quotes make a python string, inner make a literal quote in the string 
                elif 'low_pressure' in base_name:
                    call_str_matlab += "'low'"
                else:
                    call_str_matlab += "'default'"

                if stride_found:
                    call_str_matlab += (", " + stride)

                call_str_matlab += ')"'

                print("call_str_matlab = ", call_str_matlab)

                code = subprocess.call(call_str_matlab, shell=True)
                if code is None:
                    print('call returned prematurely')

                left = base_name
                right = 'flow_and_pressure_'
                base_name_new = base_name + '_flow_trace'
                nsteps = 963

                call_str_concatenate = 'python ~/copies_scripts/concatenate_jpegs.py '
                call_str_concatenate += left + ' ' + right + ' ' + base_name_new + ' ' + str(nsteps)

                code = subprocess.call(call_str_concatenate, shell=True)
                if code is None:
                    print('call returned prematurely')


    # ffmpeg from here so variables are all in place
    movie_string = 'ffmpeg -y -framerate 60 -i '
    movie_string += base_name_new
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name_new + '.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print('something wrong in movie make, call returned prematurely')


    # # reduce by 10x
    # # 60 input, 60 output is 10x slow motion
    # # 600 input, 60 output is real time 

    movie_string = 'ffmpeg -y -framerate 600 -i '
    movie_string += base_name_new
    movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
    movie_string += base_name_new + '_real_time.mp4'

    code = subprocess.call(movie_string, shell=True)
    if code is None:
        print('something wrong in movie make, call returned prematurely')
