
import subprocess
import os


cwd = os.getcwd()


session_file = '/scratch/adk354/three_slice_128_high_res.session'
# session_file = 'three_slice_128_high_res.session'

lagrangian_visit_file = cwd + '/lag_data.visit'
eulerian_visit_file   = cwd + '/dumps.visit'

print 'trying to open ', lagrangian_visit_file
print 'trying to open ', eulerian_visit_file

data = (lagrangian_visit_file, eulerian_visit_file)

RestoreSessionWithDifferentSources(session_file, 0, data)


# get some output names
cwd_split = cwd.split('/')

# name files after the job if easy
# path is always /home/adk354/scratch/JOB_NAME
if (cwd_split[0]   == '') and (cwd_split[1]   == 'scratch') and (cwd_split[0]   == 'adk354') and (len(cwd_split) >= 4):
    base_name = cwd_split[3]
else:
    base_name = 'frames'

s = SaveWindowAttributes()
s.fileName = base_name
s.outputDirectory = cwd
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

movie_string = 'ffmpeg -framerate 30 -i '
movie_string += base_name
movie_string += '%4d.jpeg -vf scale=iw*.25:ih*.25 -r 30 -c:v libx264 -preset veryslow -crf 18 '
movie_string += base_name + '.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

quit()
