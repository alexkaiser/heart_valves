
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
    print 'defaulting to default file name'
    base_name = "mitral_tree"
else: 
    base_name = str(sys.argv[1])


# make a mesh plot for dumb reasons
AddPlot("Mesh", base_name + "_vertices")
DrawPlots()

exp_db = ExportDBAttributes() 

exp_db.db_type = "XYZ" 

exp_db.variables = (base_name + "_mesh")


# ExportDatabase(exp_db)  


for state in range(TimeSliderGetNStates()):
    SetTimeSliderState(state)
    print 'state = ', state

    exp_db.filename = base_name + "_lines3d_" + str('%010d' % state)   
    ExportDatabase(exp_db)  
    
    #if state >= 0:
    #     break
    
    
print 'script cleared without crash'
 
quit() 
 
