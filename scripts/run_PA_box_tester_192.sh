#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=aortic_1
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden,willhies
# #SBATCH --exclude=sh-107-[59-64]

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

# executable is here 
SRCDIR=$PWD

# run in scratch, name with the job name
RUNDIR=$SCRATCH/pa_${SLURM_JOBID/.*}_192_a9aa59d_r_250_to_exp_pressure_adjusted_half_rv_rho_mu_updated_g_on
mkdir $RUNDIR

# set up run info 
BASE_NAME=aortic_no_partition_192
BASE_NAME_VESSEL=vessel
INPUT_NAME=input_PA_box_tester_valve_192
RUN_LINE="mpiexec --bind-to core -report-bindings main_rv_pa"
OPTIONS="-velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_max_it 1 -velocity_ksp_norm_type none > output.txt 2>&1"
SESSION_NAME="pa_valve_tester.session"
VIEW_CLIPPING="0.2"

pwd

# move stuff the the run directory 
cd $SRCDIR
cp left_pa_bdry.vertex                $RUNDIR
cp right_pa_bdry.vertex               $RUNDIR 
cp right_ventricle_bdry.vertex        $RUNDIR
cp fourier_coeffs*pa*                 $RUNDIR
cp fourier_coeffs_right_ventricle.txt $RUNDIR
cp $BASE_NAME*                        $RUNDIR
cp $BASE_NAME_VESSEL*                 $RUNDIR
cp $INPUT_NAME                        $RUNDIR
cp *.cpp                              $RUNDIR
cp main_rv_pa                         $RUNDIR
cp watchdog_job_restart.py            $RUNDIR
cp kill_all_mpi.sh                    $RUNDIR
cp run_PA_box_tester_192.sh           $RUNDIR
cp run_parallel_movie.py              $RUNDIR

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

# convert to paraview formats 
visit -cli -nowin -s ~/copies_scripts/run_parallel_convert_visit_to_paraview.py $SLURM_NTASKS $SLURM_NTASKS

sbatch ~/copies_scripts/post_process_pa.sh 

