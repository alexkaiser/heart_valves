
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

t = 0.0

prefix = '''{
  "file-series-version" : "1.0",
  "files" : [
'''

suffix = '''  ]
}
'''

for n in range(nsteps):
    for proc in range(nprocs):

        if not initialized:

            filename_out = basename + '.' + extension + '.series'
            print "filename_out = ", filename_out

            f_write = open(filename_out, 'w')
            f_write.write(prefix)
            initialized = True

        tmp_str  = '    { '

        tmp_str += '"name" : "'
        tmp_str += basename + str(n).zfill(4) + '.' 
        if nprocs > 1:
            tmp_str += str(proc) + '.'        
        tmp_str += extension + '", '

        tmp_str += '"time" : ' + '{:.14f}'.format(t) + ', '

        if nprocs > 1:
            tmp_str += '"part" : ' + str(proc)
        
        tmp_str += ' }'

        if not ((n == (nsteps-1)) and (proc == (nprocs-1))):
            tmp_str += ','      # trailing comma not tolerated at end of last line 

        tmp_str += '\n'

        f_write.write(tmp_str)

    t += dt  

f_write.write(suffix)
f_write.close()
