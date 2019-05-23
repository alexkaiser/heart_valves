#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=167:00:00
#SBATCH --mem=125GB
#SBATCH --job-name=mitral_256
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --exclusive


module purge
module load openmpi/gnu/2.0.2

# executable is here 
SRCDIR=$HOME/mitral_fully_discrete

# run in scratch, name with the job name
RUNDIR=$SCRATCH/mitral_cycle_PERIODIC_${SLURM_JOBID/.*}_256_58556997_mean_reverting_convergences
mkdir $RUNDIR

# set up run info 
BASE_NAME=mitral_tree_256
INPUT_NAME=input_mitral_tree_cycle_256_PERIODIC_systole_skeleton
RUN_LINE="mpiexec --bind-to core -report-bindings main3d"
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
cp run_cycle_N256_PERIODIC_systolic.sbatch   $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS"
