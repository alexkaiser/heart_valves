#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=three_layer_thicken
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

# executable is here 
SRCDIR=$HOME/mitral_fully_discrete

# run in scratch, name with the job name
RUNDIR=$SCRATCH/mitral_lv_hocm_d_${SLURM_JOBID/.*}_256_67b0ae80_three_layer_thicken_ds_pt05cm
mkdir $RUNDIR

# set up run info 
BASE_NAME=mitral_no_partition_256
BASE_NAME_LV=LV_pt05cm_hocm_d
INPUT_NAME=input_mitral_LV_256_hocm_d
RUN_LINE="mpiexec --bind-to core -report-bindings main_with_lv"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME*                               $RUNDIR
cp $BASE_NAME_LV*                            $RUNDIR
cp aorta_bdry.vertex                         $RUNDIR
cp atrium_bdry.vertex                        $RUNDIR
cp $INPUT_NAME                               $RUNDIR
cp *.cpp                                     $RUNDIR
cp main_with_lv                              $RUNDIR
cp fourier_coeffs*.txt                       $RUNDIR
cp watchdog_job_restart.py                   $RUNDIR
cp kill_all_mpi.sh                           $RUNDIR
cp run_mitral_lv_hocmd_N256_sherlock.sbatch  $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS"
