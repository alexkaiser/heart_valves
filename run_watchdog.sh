#!/bin/bash
#PBS -l nodes=2:ppn=20
#PBS -l walltime=4:00:00
#PBS -l mem=124GB
#PBS -N mitral_cycle
#PBS -M kaiser@cims.nyu.edu
#PBS -m abe
#PBS -j oe
#PBS -n

# Copyright (c) 2019, Alexander D. Kaiser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

module purge
module load openmpi/gnu/1.6.5

# executable is here 
SRCDIR=$HOME/mitral_fully_discrete

# run in scratch, name with the job name
RUNDIR=$SCRATCH/mitral_cycle_PERIODIC_${PBS_JOBID/.*}
mkdir $RUNDIR

# set up run info 
BASE_NAME=mitral_tree_128
INPUT_NAME=input_mitral_tree_cycle_128_PERIODIC
RUN_LINE="mpirun -hostfile hostfile --bind-to-core -report-bindings main3d"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME.*                 $RUNDIR
cp $INPUT_NAME                  $RUNDIR
cp main.cpp                     $RUNDIR
cp main3d                       $RUNDIR
cp fourier_coeffs.txt           $RUNDIR
cp run_watchdog.sh              $RUNDIR
cp watchdog_job_restart.py      $RUNDIR
cp kill_all_mpi.sh              $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# sets the nodes to be n per node 
export procs_per_node=20

hostfile=hostfile
{
    for node in $(cat $PBS_NODEFILE | uniq); do
	for((i=0; i<$procs_per_node; i++)); do
	    echo $node
	done
    done
} > $hostfile


# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS"




