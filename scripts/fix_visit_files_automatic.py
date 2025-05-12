import glob
from natsort import natsorted

lag_dirs = natsorted(glob.glob("lag_data.cycle_*"))
eulerian_dirs = natsorted(glob.glob("visit_dump.*"))

dumps = open('dumps.visit', 'w')
lag   = open('lag_data.visit', 'w')

for lag_dir in lag_dirs:
    lag.write(lag_dir + '/' + lag_dir + '.summary.silo\n') 

for visit_dir in eulerian_dirs:
    dumps.write(visit_dir +  '/summary.samrai\n')

dumps.close()
lag.close()

