from __future__ import print_function

print("in script")

OpenDatabase("dumps.visit")

print("open passed")

f = open('times.txt', 'w')

SetTimeSliderState(0)
Query('Time')
# if GetQueryOutputValue() != 0.0:
#     raise ValueError('First timestep does not have zero time')
# else: 
#     print("found zero time" )

for state in range(TimeSliderGetNStates()): 

    SetTimeSliderState(state)
    Query('Time')
    t  = GetQueryOutputValue()
    print("t = " , t)
    f.write(  "{:0.10e}".format(t) + "\n")

f.close()
