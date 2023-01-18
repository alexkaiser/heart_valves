from __future__ import print_function
import subprocess
import os
import sys 

if len(sys.argv) < 2:
    raise InputError("Must provide a script name")

session_name = sys.argv[1]

if len(sys.argv) >= 4:
    n_procs  = int(sys.argv[2])
    proc_num = int(sys.argv[3])
else:
    n_procs = 1
    proc_num = 0

if len(sys.argv) >= 5:
    try:
        view_clipping = float(sys.argv[4])
    except ValueError:
        print ("Not a float")
        view_clipping = None 
else:
    view_clipping = None

cwd = os.getcwd()

session_file = session_name
#'~/mitral_fully_discrete/normal_1_two_slice_two_valve.session'
# session_file = 'three_slice_128_high_res.session'

lagrangian_visit_file = cwd + '/lag_data.visit'
eulerian_visit_file   = cwd + '/dumps.visit'

print ('trying to open ', lagrangian_visit_file)
print ('trying to open ', eulerian_visit_file)

# data = (lagrangian_visit_file, eulerian_visit_file)
data = (eulerian_visit_file, lagrangian_visit_file)

RestoreSessionWithDifferentSources(session_file, 0, data)
# RestoreSession(session_file, 0)

print ('restore passed')

# get some output names
cwd_split = cwd.split('/')

# name files after the job if easy
# path is always /home/adk354/scratch/JOB_NAME
if (len(cwd_split) >= 7) and (cwd_split[1] == "expanse"):
    base_name = cwd_split[6]
elif (len(cwd_split) >= 5):
    base_name = cwd_split[4]
else:
    base_name = 'frames'

s = GetSaveWindowAttributes()
s.fileName = base_name
s.outputDirectory = cwd
s.format = s.JPEG
s.family = 0
#s.saveTiled = 1
#s.width = 1920*4
SetSaveWindowAttributes(s)


for state in range(TimeSliderGetNStates()):
    # quick hack parallelism 
    if (state % n_procs) == proc_num:
        try:
            if state % (TimeSliderGetNStates()/20) == 0:
                print ('On state ', state, 'of ', TimeSliderGetNStates())
        except: 
            # ignore message if it gives a divide by zero 
            # since there are not enough states 
            pass

        SetTimeSliderState(state)

        s.fileName = base_name + str(state).zfill(4) + ".jpeg"
        SetSaveWindowAttributes(s)

        if view_clipping is not None:
            view_obj = GetView3D()
            view_obj.nearPlane = view_clipping
            SetView3D(view_obj)

        # make sure to update in time
        DrawPlots()

        if view_clipping is not None:
            view_obj = GetView3D()
            view_obj.nearPlane = view_clipping
            SetView3D(view_obj)

        SaveWindow()


# if proc_num == 0:

#     # # call ffmpeg from here so variables are all in place
#     # code = subprocess.call('module load ffmpeg', shell=True)
#     # if code is None:
#     #     print 'ffmpeg module load did not finish'
#     #     print 'weird results likely'

#     movie_string = 'ffmpeg -framerate 60 -i '
#     movie_string += base_name
#     movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
#     movie_string += base_name + '.mp4'

#     code = subprocess.call(movie_string, shell=True)
#     if code is None:
#         print 'something wrong in movie make, call returned prematurely'


#     # # reduce by 10x
#     # # 60 input, 60 output is 10x slow motion
#     # # 600 input, 60 output is real time 

#     # movie_string = 'ffmpeg -framerate 600 -i '
#     # movie_string += base_name
#     # movie_string += '%4d.jpeg -vf scale=1920:-2 -r 60 -c:v libx264 -preset veryslow -crf 18 '
#     # movie_string += base_name + '_real_time.mp4'

#     # code = subprocess.call(movie_string, shell=True)
#     # if code is None:
#     #     print 'something wrong in movie make, call returned prematurely'

quit()
