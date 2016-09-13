
import sys 

if len(sys.argv) < 3:
    print 'Cannot update without stride and final data number'
else: 
    stride   = int(sys.argv[1])
    max_step = int(sys.argv[2])

    dumps = open('dumps.visit', 'w')
    lag   = open('lag_data.visit', 'w')

    state = 0 
    dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
    lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 


    for state in range(stride, max_step, stride):
        dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
        lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 

    state = max_step
    dumps.write('visit_dump.' + str('%05d' % state)    +  '/summary.samrai\n')
    lag.write('lag_data.cycle_' + str('%06d' % state) + '/lag_data.cycle_' + str('%06d' % state) + '.summary.silo\n') 

    dumps.close()
    lag.close()

