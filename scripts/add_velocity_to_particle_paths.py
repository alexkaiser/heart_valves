import os
import vtk
import sys
import numpy as np
from vtk.util import numpy_support

def read_vtu(infile):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(infile)
    reader.Update()
    data = reader.GetOutput()
    return data

def write_vtu(outfile,outdata):
    #Write the vtp model to file
    modelwrite=vtk.vtkXMLUnstructuredGridWriter()
    #modelwrite.SetInputData(model)
    modelwrite.SetInputData(outdata)
    print("Writing appended file..." + outfile)
    modelwrite.SetFileName(outfile)
    modelwrite.Write()
    return

if len(sys.argv) < 5:
    raise TypeError('required agrs: basename nsteps dt') 
else: 
    basename = sys.argv[1]
    extension = sys.argv[2]
    nsteps = int(sys.argv[3])
    dt = float(sys.argv[4])

if extension != 'vtu':
    raise TypeError('unsupported extension')

#Set file suffixes
first_timestep=0

#Set velocity to 0 at all points for first timestep
FileName=basename+("%04d"%first_timestep)+".vtu"
data1 = read_vtu(FileName)
coords1 = data1.GetPoints().GetData()

velocity = numpy_support.numpy_to_vtk(np.zeros(np.array(coords1).shape), deep=True, array_type=vtk.VTK_DOUBLE)
velocity.SetName("velocity")
data1.GetPointData().AddArray(velocity)

write_vtu(basename + "_vel"+("%04d"%first_timestep)+".vtu", data1)

for i in range(1,nsteps):

    FileName=basename+("%04d"%(i-1))+".vtu"
    data1 = read_vtu(FileName)
    coords1 = data1.GetPoints().GetData()
    FileName=basename+("%04d"%(i  ))+".vtu"
    data2 = read_vtu(FileName)
    coords2 = data2.GetPoints().GetData()
    
    velocity = numpy_support.numpy_to_vtk((1.0/dt) * (np.array(coords2)-np.array(coords1)), deep=True, array_type=vtk.VTK_DOUBLE)
    velocity.SetName("velocity")
    data2.GetPointData().AddArray(velocity)

    write_vtu(basename + "_vel"+("%04d"%(i))+".vtu", data2)

