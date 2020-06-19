
import sys 

if len(sys.argv) < 5:
    raise TypeError('required agrs: basename extension nsteps dt') 
else: 
    basename = sys.argv[1]
    extension = sys.argv[2]
    nsteps = int(sys.argv[3])
    dt = float(sys.argv[4])

if len(sys.argv) >= 6:
    nprocs = int(sys.argv[5])
else: 
    nprocs = 1 

initialized = False 

# example 
'''
<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1"
         byte_order="LittleEndian"
         compressor="vtkZLibDataCompressor">
  <Collection>
    <DataSet timestep="0" group="" part="0"
             file="examplePVD/examplePVD_T0000.vtp"/>
    <DataSet timestep="1" group="" part="0"
             file="examplePVD/examplePVD_T0001.vtp"/>
    <DataSet timestep="2" group="" part="0"
             file="examplePVD/examplePVD_T0002.vtp"/>
    <DataSet timestep="3" group="" part="0"
             file="examplePVD/examplePVD_T0003.vtp"/>
    <DataSet timestep="4" group="" part="0"
             file="examplePVD/examplePVD_T0004.vtp"/>
  </Collection>
</VTKFile>
'''




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
