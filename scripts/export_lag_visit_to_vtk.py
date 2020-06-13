
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

print 'we out here'

OpenDatabase("lag_data.visit")
 
print 'data base open line passed'

if len(sys.argv) <= 1:
    raise TypeError('required agrs: base_name')
else: 
    base_name = sys.argv[1]

if len(sys.argv) >= 3:
    extension = sys.argv[2]
else:
    print "using default extension vtu"
    extension = 'vtu'

if len(sys.argv) >= 6:
    nprocs = int(sys.argv[4])
    proc_num = int(sys.argv[5])
else: 
    print "using default proc_num 0, nprocs = 1"
    proc_num = 0
    nprocs = 1

# make a mesh plot for dumb reasons
AddPlot("Mesh", base_name)
DrawPlots()

exp_db = ExportDBAttributes() 
exp_db.db_type = 'VTK' 
exp_db.variables = (base_name)

export_opts = GetExportOptions("VTK")
print "export_opts = ", export_opts

export_opts['Binary format'] = 1
export_opts['XML format'] = 1


# also write a .series file for paraview 
initialized = False 

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

    if not initialized:

        filename_out = base_name + '.' + extension + '.series'
        print "filename_out = ", filename_out

        f_write = open(filename_out, 'w')
        f_write.write(prefix)
        initialized = True

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
print 'script cleared without crash'
 
quit() 
