#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=aortic_1
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden
#SBATCH --exclude=sh03-16n02
# #SBATCH --exclude=sh-107-[59-64]

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

TOTAL_TASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))

# executable is here 
SRCDIR=$PWD

# run in scratch, name with the job name
RUNDIR=$SCRATCH/aortic_${SLURM_JOBID/.*}_192_8b6ec8_initial_quad_c30_r20_loaded
mkdir $RUNDIR

# set up run info 
BASE_NAME=aortic_no_partition_192
BASE_NAME_VESSEL=aorta_192
INPUT_NAME=input_aortic_192_graft
RUN_LINE="srun main_aorta"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"
SESSION_NAME="aortic_384_with_aorta.session"
VIEW_CLIPPING="-0.2"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp $BASE_NAME*                               $RUNDIR
cp $BASE_NAME_VESSEL*                        $RUNDIR
cp cylinder_bdry_max_192.vertex              $RUNDIR
cp cylinder_bdry_min_192.vertex              $RUNDIR
cp $INPUT_NAME                               $RUNDIR
cp *.cpp                                     $RUNDIR
cp main_aorta                                $RUNDIR
cp fourier_coeffs_ventricle.txt              $RUNDIR
cp fourier_coeffs_Q_mi.txt                   $RUNDIR
cp fourier_coeffs_lv_activation_two_hill.txt $RUNDIR
cp watchdog_job_restart.py                   $RUNDIR
cp kill_all_mpi.sh                           $RUNDIR
cp run_aortic_192_graft.sbatch               $RUNDIR
cp run_parallel_movie.py                     $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
env_log=$RUNDIR/env.log
rm -rf $env_log
env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log

# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS" 

# touch done.txt

# load stuff for movie making 
source ~/.bash_profile
visit -cli -nowin -s ~/heart_valves/scripts/run_parallel_convert_visit_to_paraview.py $TOTAL_TASKS $TOTAL_TASKS

# sbatch ~/heart_valves/scripts/post_process_aorta.sh 
source ~/.bash_profile

sbatch ~/copies_scripts/post_process_graft.sh


