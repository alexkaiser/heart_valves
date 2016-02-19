
# python visit script

# looks for a file named "lag_data.visit"
# must be in the current directory 
# exports this to a simple three coordinate list for all times 

# to run 

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
 
