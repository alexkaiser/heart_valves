
import subprocess
import os

# Copyright (c) 2019, Alexander D. Kaiser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
