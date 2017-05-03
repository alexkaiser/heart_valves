
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

if (n != 128) and (n != 256) and (n != 512):
    print 'session files only available for 128,256,512 resolutions so far'
    quit()


session_file = '/scratch/adk354/three_slice_' + str(n) + '_high_res.session'
# session_file = 'three_slice_128_high_res.session'

lagrangian_visit_file = cwd + '/lag_data.visit'
eulerian_visit_file   = cwd + '/dumps.visit'

print 'trying to open ', lagrangian_visit_file
print 'trying to open ', eulerian_visit_file

data = (lagrangian_visit_file, eulerian_visit_file)

RestoreSessionWithDifferentSources(session_file, 0, data)

print 'restore passed'

# get some output names
cwd_split = cwd.split('/')

# name files after the job if easy
# path is always /home/adk354/scratch/JOB_NAME
if (cwd_split[0]   == '') and (cwd_split[1]   == 'scratch') and (cwd_split[2]   == 'adk354') and (len(cwd_split) >= 4):
    base_name = cwd_split[3]
else:
    base_name = 'frames'

s = SaveWindowAttributes()
s.fileName = base_name
s.outputDirectory = cwd
s.format = s.JPEG
s.saveTiled = 1
s.width = 1920*4
SetSaveWindowAttributes(s)


for state in range(TimeSliderGetNStates()):
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

movie_string = 'ffmpeg -framerate 60 -i '
movie_string += base_name
movie_string += '%4d.jpeg -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 18 '
movie_string += base_name + '.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'


# reduce by 10x
# 60 input, 60 output is 10x slow motion
# 600 input, 60 output is real time 

movie_string = 'ffmpeg -framerate 600 -i '
movie_string += base_name
movie_string += '%4d.jpeg -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 18 '
movie_string += base_name + '_real_time.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

quit()
