
# python visit script

# looks for a file named "lag_data.visit"
# must be in the current directory 
# exports this to a simple three coordinate list for all times 

# to run 


print 'we out here'

OpenDatabase("lag_data.visit")
 
print 'data base open line passed'


# make a mesh plot for dumb reasons
AddPlot("Mesh", "mitral_tree_vertices")
DrawPlots()

exp_db = ExportDBAttributes() 

exp_db.db_type = "XYZ" 

exp_db.variables = ("mitral_tree_mesh")


# ExportDatabase(exp_db)  


for state in range(TimeSliderGetNStates()):
    SetTimeSliderState(state)
    print 'state = ', state

    exp_db.filename = "mitral_mesh_lines3d_" + str('%010d' % state)   
    ExportDatabase(exp_db)  
    
    if state >= 0:
         break
    
    
print 'script cleared without crash'
 
 
 
