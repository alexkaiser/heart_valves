
# python visit script

# looks for a file named "dumps.visit"
# must be in the current directory 
# 

# import matplotlib 


print 'we out here'



OpenDatabase("dumps.visit")
 
print 'data base open line passed'

# just make sure it's cleaned up 
DeleteAllPlots()




AddPlot("Pseudocolor", "U_z")
AddOperator("Slice")
# SetActivePlots(plot)

SetActivePlots(1)

# print "supported queries: ", Queries()


slice = SliceAttributes()

# print 'slice attributes, before modification'
# print slice
# print 'slice.originType = ', slice.originType, 'its type = ', type(slice.originType)

slice.originType  = 0  # point, required to be an int
slice.originPoint = (0, 0, 9) 
slice.project2d   = 0
slice.normal      = (0, 0, -1)
# print 'after setting, but before call to set options'
# print slice
SetOperatorOptions(slice)
# print 'after setting and call to set options'
# print slice

DrawPlots()

# exp_db = ExportDBAttributes() 
# exp_db.db_type = "XYZ" 
# exp_db.variables = ("mitral_tree_mesh")

# ExportDatabase(exp_db)  

Query('Variable Sum')
sum = GetQueryOutputValue()
print 'sum = ', sum  




# geometry parameters 
output_stride = 10
dt = output_stride * 0.00002

# x,y domain is [-L, L] total 
L = 3.0
N = 32
dx = 2.0 * L / N 

times    = []
flux     = []
net_flux = []

t          = 0.0
flux_prev  = 0.0 
total_flux = 0.0


for state in range(TimeSliderGetNStates()):
    
    if state % (TimeSliderGetNStates()/20) == 0:
        print 'On state ', state, 'of ', TimeSliderGetNStates()
    
    SetTimeSliderState(state)
    times.append(t)

    DrawPlots()
    
    # add the z velocity components in the slice 
    Query('Variable Sum')
    
    # multiplying by dx^2 gives a second order midpoint approx 
    # to the current flux
    # add a sign because we want inward flux  
    flux_current = -dx * dx * GetQueryOutputValue()

    flux.append(flux_current)

    # simple trapezoidal rule for net flux 
    total_flux += 0.5 * dt * (flux_current + flux_prev)
    net_flux.append(total_flux)
    
    # variable update for next step 
    t += dt 
    flux_prev = flux_current
    
    #if state > 200:
    #     break

    


f = open('plot_flux.m', 'w')

f.write('times = ')
f.write( str(times) )
f.write(';\n\n')

f.write('flux = ')
f.write( str(flux) )
f.write(';\n\n')

f.write('net_flux = ')
f.write( str(net_flux) )
f.write(';\n\n')

plot_code = '''
fig = figure; 
plot(times, flux); 
title('flux over time'); 
xlabel('t'); 
ylabel('flux (cm^3 / s)  '); 

fig = figure; 
plot(times, net_flux); 
title('nex flux over time'); 
xlabel('t'); 
ylabel('net flux (cm^3)'); 
'''

f.write(plot_code) 




f.close()

    
print 'script cleared without crash'
 
 
 


'''
# set the active plot to the mesh plot 
# SetActivePlots(1)
AddPlot("Mesh", "amr_mesh")

options = GetPlotOptions()
print 'plot options = ', options

mesh_attributes = GetMeshManagementAttributes()
print "mesh_attributes = ", mesh_attributes
print "discretizationTolernace = ", mesh_attributes.discretizationTolerance
'''
'''
box = BoxAttributes()
print 'box = ', box 
''' 
