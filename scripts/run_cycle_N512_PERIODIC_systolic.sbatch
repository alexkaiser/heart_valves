#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=48:00:00
#SBATCH --mem=250GB
#SBATCH --job-name=mitral_cycle
#SBATCH --mail-user=kaiser@cims.nyu.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --exclude=sh03-16n02

module purge
module load openmpi/gnu/2.0.2

# executable is here 
SRCDIR=$HOME/heart_valves

# run in scratch, name with the job name
RUNDIR=$SCRATCH/mitral_cycle_PERIODIC_${SLURM_JOBID/.*}_512_git_9f010d40_best_general_setup_at_high_res
mkdir $RUNDIR

# set up run info 
BASE_NAME=mitral_tree_512
INPUT_NAME=input_mitral_tree_cycle_512_PERIODIC_systole_skeleton
RUN_LINE="srun main3d"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME.*                              $RUNDIR
cp $INPUT_NAME                               $RUNDIR
cp *.cpp                                     $RUNDIR
cp main3d                                    $RUNDIR
cp fourier_coeffs*.txt                       $RUNDIR
cp watchdog_job_restart.py                   $RUNDIR
cp kill_all_mpi.sh                           $RUNDIR
cp run_cycle_N512_PERIODIC_systolic.sbatch   $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS"
