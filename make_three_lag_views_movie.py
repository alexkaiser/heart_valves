
import subprocess
import os


cwd = os.getcwd()

# find out what power of two we are working with
count = 1
n = 2 
while True: 
    if os.path.isfile('../mitral_tree_' + str(n) + '.vertex'): 
        base = n
        break

    n *= 2
    count += 1
    if count > 20:
        break

if (n != 128) and (n != 256):
    print 'session files only available for 128,256 resolutions so far'


session_file = '/scratch/adk354/three_views_lagrangian_' + str(n) + '.session'
# session_file = 'three_slice_128_high_res.session'

lagrangian_visit_file = cwd + '/lag_data.visit'

print 'trying to open ', lagrangian_visit_file

data = (lagrangian_visit_file)

RestoreSessionWithDifferentSources(session_file, 0, data)

print 'restore passed'

# get some output names
cwd_split = cwd.split('/')

# name files after the job if easy
# path is always /home/adk354/scratch/JOB_NAME
if (cwd_split[0]   == '') and (cwd_split[1]   == 'scratch') and (cwd_split[2]   == 'adk354') and (len(cwd_split) >= 4):
    base_name = cwd_split[3] + '_lagrangian'
else:
    base_name = 'frames' + '_lagrangian'

s = SaveWindowAttributes()
s.fileName = base_name
s.outputDirectory = cwd
s.saveTiled = 1
s.width = 1920*4
SetSaveWindowAttributes(s)


for state in range(0,TimeSliderGetNStates(),2):
    try:
        if state % (TimeSliderGetNStates()/20) == 0:
            print 'On state ', state, 'of ', TimeSliderGetNStates()
    except: 
        # ignore message if it gives a divide by zero 
        # since there are not enough states 
        pass

    SetTimeSliderState(state)

    # make sure to update in time
    DrawPlots()

    SaveWindow()

# call ffmpeg from here so variables are all in place
code = subprocess.call('module load ffmpeg', shell=True)
if code is None:
    print 'ffmpeg module load did not finish'
    print 'weird results likely'

# took output stride of two
movie_string = 'ffmpeg -framerate 60 -i '
movie_string += base_name
movie_string += '%4d.jpeg -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 18 '
movie_string += base_name + '.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'


code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

quit()