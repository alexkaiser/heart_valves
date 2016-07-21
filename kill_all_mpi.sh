#!/bin/bash

# read in the old environment
env_log=env.log
while read -r line; do
    export "$line"
done < $env_log

{
for node in $(cat $PBS_NODEFILE | uniq); do
    ssh -x $node "killall -9 main3d; killall -9 mpirun; killall -9 orted"
done
}
