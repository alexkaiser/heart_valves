
import subprocess
import os

# get some output names
cwd = os.getcwd()
cwd_split = cwd.split('/')

# name files after the job if easy
# path is always /home/adk354/scratch/JOB_NAME
if (cwd_split[0]   == '') and (cwd_split[1]   == 'scratch') and (cwd_split[2]   == 'adk354') and (len(cwd_split) >= 4):
    base_name = cwd_split[3]
else:
    base_name = 'frames'




# call ffmpeg from here so variables are all in place
code = subprocess.call('module load ffmpeg', shell=True)
if code is None:
    print 'ffmpeg module load did not finish'
    print 'weird results likely'

movie_string = 'ffmpeg -framerate 60 -i particle_views_%6d_panels.jpg'
movie_string += ' -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 35.3 '
movie_string += base_name + '_particles.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

movie_string = 'ffmpeg -framerate 60 -i particle_views_%6d_panels.jpg'
movie_string += ' -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 17 '
movie_string += base_name + '_uncompressed_particles.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'


# movie_string = 'ffmpeg -framerate 60 -i particle_views_%6d_panels.jpg'
# movie_string += ' -vf scale=iw*.25:ih*.25 -r 30 -c:v libx264 -preset veryslow -crf 33 '
# movie_string += base_name + '_particles_30fps.mp4'

# code = subprocess.call(movie_string, shell=True)
# if code is None:
#     print 'something wrong in movie make, call returned prematurely'

# reduce by 10x
# 60 input, 60 output is 10x slow motion
# 600 input, 60 output is real time 

movie_string = 'ffmpeg -framerate 600 -i particle_views_%6d_panels.jpg'
movie_string += ' -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 24 '
movie_string += base_name + '_particles_real_time.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

movie_string = 'ffmpeg -framerate 600 -i particle_views_%6d_panels.jpg'
movie_string += ' -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 17 '
movie_string += base_name + '_uncompressed_particles_real_time.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

# movie_string = 'ffmpeg -framerate 600 -i particle_views_%6d_panels.jpg'
# movie_string += ' -vf scale=iw*.25:ih*.25 -r 30 -c:v libx264 -preset veryslow -crf 21 '
# movie_string += base_name + '_particles_real_time_30fps.mp4'

# code = subprocess.call(movie_string, shell=True)
# if code is None:
#     print 'something wrong in movie make, call returned prematurely'

quit()
