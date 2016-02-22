
# python visit script

# looks for a file named "dumps.visit"
# must be in the current directory 
# 


print 'we out here'

OpenDatabase("dumps.visit")
 
print 'data base open line passed'

# just make sure it's cleaned up 
DeleteAllPlots()


AddPlot("Pseudocolor", "U_z")
AddOperator("Slice")
SetActivePlots(1)

# print "supported queries: ", Queries()


slice = SliceAttributes()
# print 'slice attributes, before modification'
# print slice
# print 'slice.originType = ', slice.originType, 'its type = ', type(slice.originType)
slice.originType  = 0  # point, required to be an int
slice.originPoint = (0, 0, 3) 
slice.project2d   = 0
slice.normal      = (0, 0, 1)
SetOperatorOptions(slice)


DrawPlots()


# Using this is bad, sometimes visit considers different 
# numbers of cells to be "in" the slice, 
# depending on alignment in some complicated way
# Use average instead 
# Query('Variable Sum')

# get the dimensions 
Query("SpatialExtents", use_actual_data=0)
coords = GetQueryOutputValue()
x_width = coords[1] - coords[0]
y_width = coords[3] - coords[2]

print 'x_width = ', x_width 
print 'y_width = ', y_width 

# x,y domain is [-L, L] total 
len_total = x_width

# make sure things are equal, otherwise need to fix 
assert y_width == len_total
    
area = len_total**2
print 'area = ', area 

Query("NumNodes", use_actual_data=1)
total_nodes_in_slice = GetQueryOutputValue()
print 'total_nodes_in_slice = ', total_nodes_in_slice


times    = []
flux     = []
net_flux = []

Query('Time')
t          = GetQueryOutputValue()
t_prev     = t 
flux_prev  = 0.0 
total_flux = 0.0


for state in range(TimeSliderGetNStates()):
    
    try: 
        if state % (TimeSliderGetNStates()/20) == 0:
            print 'On state ', state, 'of ', TimeSliderGetNStates()
    except: 
        # ignore message if it gives a divide by zero 
        # since there are not enough states 
        pass 
    
    t_prev = t 
    Query('Time')
    t  = GetQueryOutputValue()
    dt = t - t_prev
        
    SetTimeSliderState(state)
    times.append(t)

    # Don't need to draw the plot to get the flux 
    # DrawPlots()
    
    # add the z velocity components in the slice 
    Query('Average Value')
    
    # Average value gives second order midpoint 
    #     to the current average value
    # Multiply by area to get integral 
    # Add a sign because we want inward flux  
    avg = GetQueryOutputValue()
    flux_current = -area * avg  

    flux.append(flux_current)
    print 'flux = ', flux_current 

    # simple trapezoidal rule for net flux 
    total_flux += 0.5 * dt * (flux_current + flux_prev)
    net_flux.append(total_flux)
    
    # variable update for next step 
    flux_prev = flux_current
    


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
 
 
 

