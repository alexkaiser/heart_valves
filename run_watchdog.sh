#!/bin/bash
#PBS -l nodes=3:ppn=20
#PBS -l walltime=4:00:00
#PBS -l mem=186GB
#PBS -N mitral_cycle
#PBS -M kaiser@cims.nyu.edu
#PBS -m abe
#PBS -j oe
#PBS -n

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
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output_20_per_node.txt 2>&1"

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
n=20
hostfile=hostfile
{
    for node in $(cat $PBS_NODEFILE | uniq); do
	for((i=0; i<$n; i++)); do
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




