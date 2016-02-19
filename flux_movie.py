
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
SetActivePlots(1)


p = PseudocolorAttributes()

print p 

p.min     = -10.0
p.minFlag = 1 
p.max     = 10.0
p.maxFlag = 1
SetPlotOptions(p)




# slice and setup 
AddOperator("Slice")
slice = SliceAttributes()
slice.originType  = 0  # point, required to be an int
slice.originPoint = (0, 0, 3) 
slice.project2d   = 2
slice.normal      = (0, 0, 1)
SetOperatorOptions(slice)


DrawPlots()

Query('Variable Sum')
sum = GetQueryOutputValue()
print 'sum = ', sum  


# set things to be fairly high resolution 
s = SaveWindowAttributes()
print 'window attributes: ', s
s.width = 4096
SetSaveWindowAttributes(s)



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
    
    try: 
        if state % (TimeSliderGetNStates()/20) == 0:
            print 'On state ', state, 'of ', TimeSliderGetNStates()
    except: 
        # ignore message if it gives a divide by zero 
        # since there are not enough states 
        pass 
        
        
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

    SaveWindow()    


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

# os.system('ffmpeg -framerate 2 -i visit%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "scale=4096:trunc(ow/a/2)*2" visit_movie.mp4') 

print 'script cleared without crash'
 
 
 
 