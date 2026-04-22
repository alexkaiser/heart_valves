#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=slices_paper
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden


module purge
module load gcc/8.1.0
module load openmpi/2.0.2

source ~/.bash_profile

module load python/3.9.0

# change dir even if not needed 
cd viz_IB3d_tree_cycle_256 

TOTAL_TASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))


python3 ~/heart_valves/scripts/convert_vtu_vertices_to_csv.py

module load matlab/R2017b
cp ../aortic_*final_data.mat . 
matlab -nodesktop -nodisplay -r 'addpath ~/valve_generator; generate_cells_file; "matlab generate cells complete"; exit;'

python3 ~/heart_valves/scripts/convert_csv_with_cells.py

python3 ~/heart_valves/scripts/fix_pvd_files.py

# renders movie 
# assumes located in main sim directory 
cd .. 
SESSION_NAME_PARAVIEW="~/heart_valves/scripts/zebrafish_3_panel_movie.py"
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS

# true flag for vertical 
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS 1 
