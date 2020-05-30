
# python visit script

# looks for a file named "lag_data.visit"
# must be in the current directory 
# exports this to a simple three coordinate list for all times 

# to run 

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


import sys 
import subprocess 

print 'we out here'

OpenDatabase("dumps.visit")
 
print 'data base open line passed'

if len(sys.argv) <= 1:
    raise TypeError('required agrs: base_name')
else: 
    base_name = sys.argv[1]

if 
    nprocs = 

if len(sys.argv) >= 3:
    extension = sys.argv[2]


# make a mesh plot for dumb reasons
AddPlot("Pseudocolor", 'P')
DrawPlots()

exp_db = ExportDBAttributes() 
exp_db.db_type = 'VTK' 
exp_db.variables = ('P', 'U', 'Omega')

print "exp_db = ", exp_db
print "exp_db.opts = ", exp_db.opts
print "exp_db.opts.types = ", exp_db.opts.types

export_opts = GetExportOptions("VTK")
print "export_opts = ", export_opts

export_opts['Binary format'] = 1
export_opts['XML format'] = 1

extension = 'vtk'

# also write a .series file for paraview 
series_file_initialized = False 

prefix = '''{
  "file-series-version" : "1.0",
  "files" : [
'''

suffix = '''  ]
}
'''

nsteps = 3 # TimeSliderGetNStates()

for state in range(nsteps):
    SetTimeSliderState(state)
    print 'state = ', state

    exp_db.filename = base_name + str(state).zfill(4)
    ExportDatabase(exp_db, export_opts)  

    Query('Time')
    t  = GetQueryOutputValue()

    if state == 1:
        dt = t 

    # sort into directories, one per timestep
    subprocess.call('mkdir ' + exp_db.filename, shell=True)
    subprocess.call('mv ' + exp_db.filename + '.* ' + exp_db.filename, shell=True)

    if not series_file_initialized:

        filename_out = base_name + '.' + extension + '.series'
        print "filename_out = ", filename_out

        f_write = open(filename_out, 'w')
        f_write.write(prefix)
        series_file_initialized = True

    tmp_str  = '    { "name" : "'
    tmp_str += base_name + str(state).zfill(4) + '.' + extension
    tmp_str += '", "time" : '
    tmp_str += '{:.14f}'.format(t)
    tmp_str += ' }'
    if state != (nsteps-1):
        tmp_str += ','      # trailing comma not tolerated at end of last line 
    tmp_str += '\n'

    f_write.write(tmp_str)

f_write.write(suffix)
f_write.close()



t = 0.0

prefix = '''<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1"
         byte_order="LittleEndian"
         compressor="vtkZLibDataCompressor">
  <Collection>
'''

suffix = '''  </Collection>
</VTKFile>
'''

for n in range(nsteps):
    for proc in range(nprocs):

        if not initialized:

            filename_out = basename + '.pvd'
            print "filename_out = ", filename_out

            f_write = open(filename_out, 'w')
            f_write.write(prefix)
            initialized = True

        tmp_str = '    <DataSet timestep="'
        tmp_str += '{:.14f}'.format(t)
        tmp_str += '" group="" part="'
        tmp_str += str(proc) + '"'

        tmp_str += ' file="'
        if nprocs > 1:
            tmp_str += basename + str(n).zfill(4) + '/' # sorted into directories 
        tmp_str += basename + str(n).zfill(4) + '.' 
        if nprocs > 1:
            tmp_str += str(proc) + '.'        
        tmp_str += extension
        tmp_str += '"/>\n'

        f_write.write(tmp_str)

    t += dt  

f_write.write(suffix)
f_write.close()






print 'script cleared without crash'
 
quit() 
