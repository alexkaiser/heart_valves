#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=47:59:00
#SBATCH --mem=185GB
#SBATCH --job-name=aortic_1
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden
# #SBATCH --exclude=sh-107-[59-64]

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

# executable is here 
SRCDIR=$PWD

# run in scratch, name with the job name
RUNDIR=$SCRATCH/aortic_${SLURM_JOBID/.*}_384_36a46e_high_res_same_as_9134a0c_rad_pt24_extra_bend
mkdir $RUNDIR

# set up run info 
BASE_NAME=aortic_384
INPUT_NAME=input_aortic_384
RUN_LINE="mpiexec --bind-to core -report-bindings main3d"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"
SESSION_NAME="aortic_384_visit2pt9.session"
VIEW_CLIPPING="-0.2"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME*                               $RUNDIR
cp aorta_bdry.vertex                         $RUNDIR
cp atrium_bdry.vertex                        $RUNDIR
cp $INPUT_NAME                               $RUNDIR
cp *.cpp                                     $RUNDIR
cp main3d                                    $RUNDIR
cp fourier_coeffs*.txt                       $RUNDIR
cp watchdog_job_restart.py                   $RUNDIR
cp kill_all_mpi.sh                           $RUNDIR
cp run_aortic_N384.sbatch                    $RUNDIR
cp run_parallel_movie.py                     $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log


# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS" 

# load stuff for movie making 
source ~/.bash_profile
python run_parallel_movie.py $SESSION_NAME $SLURM_NTASKS $VIEW_CLIPPING
