#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --time=48:00:00
#SBATCH --mem=249325M
#SBATCH --job-name=aortic_1
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=compute
#SBATCH --account=csd481

source ~/.bash_profile

TOTAL_TASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))

# executable is here 
SRCDIR=$PWD

# run in scratch, name with the job name
RUNDIR=$SCRATCH/pa_${SLURM_JOBID/.*}_384_2f261b2_variable_r_lpa_20_res_5_16th
mkdir $RUNDIR

# set up run info 
BASE_NAME=aortic_no_partition_384
BASE_NAME_VESSEL=vessel_384
INPUT_NAME=input_PA_box_tester_valve_384
RUN_LINE="mpiexec --bind-to core -report-bindings main_rv_pa"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"
SESSION_NAME="pa_valve_tester_384.session"
VIEW_CLIPPING="0.2"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp left_pa_bdry.vertex                $RUNDIR
cp right_pa_bdry.vertex               $RUNDIR 
cp right_ventricle_*vertex            $RUNDIR
cp fourier_coeffs*pa*                 $RUNDIR
cp fourier_coeffs_right_ventricle*    $RUNDIR
cp fourier_coeffs_zero*               $RUNDIR
cp $BASE_NAME*                        $RUNDIR
cp $BASE_NAME_VESSEL*                 $RUNDIR
cp $INPUT_NAME                        $RUNDIR
cp *.cpp                              $RUNDIR
cp main_rv_pa                         $RUNDIR
cp watchdog_job_restart.py            $RUNDIR
cp kill_all_mpi.sh                    $RUNDIR
cp run_PA_box_tester_384_expanse.sh   $RUNDIR
cp run_parallel_movie.py              $RUNDIR

# go to the run directory before running anything 
cd $RUNDIR

# dump current environment to file 
# env_log=$RUNDIR/env.log
# rm -rf $env_log
# env | grep -v '{' | grep -v '}' | grep -v '()' | grep -v _= > $env_log

# call python script which controls and watches the run 
python watchdog_job_restart.py "$RUN_LINE" "$INPUT_NAME" "$OPTIONS" 

# load stuff for movie making 
# source ~/.bash_profile
# python run_parallel_movie.py $SESSION_NAME $TOTAL_TASKS $VIEW_CLIPPING

# # convert to paraview formats 
# visit -cli -nowin -s ~/copies_scripts/run_parallel_convert_visit_to_paraview.py $TOTAL_TASKS $TOTAL_TASKS

# if test -f done.txt; then
#     sbatch ~/copies_scripts/post_process_pa.sh
# fi
