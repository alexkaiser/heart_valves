
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


# 10x faster
movie_string = 'ffmpeg -framerate 600 -i '
movie_string += base_name
movie_string += '%4d.jpeg -vf scale=iw*.25:ih*.25 -r 60 -c:v libx264 -preset veryslow -crf 18 '
movie_string += base_name + '_real_time.mp4'

code = subprocess.call(movie_string, shell=True)
if code is None:
    print 'something wrong in movie make, call returned prematurely'

quit()
