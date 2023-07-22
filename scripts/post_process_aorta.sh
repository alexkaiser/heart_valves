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

cp ~/mitral_fully_discrete/2_aorta_remeshed_pt5mm_capped.vtp . 
cp ~/mitral_fully_discrete/4_aorta_remeshed_pt25mm_3cm_extender_layers_constriction.vtu . 
cp ~/mitral_fully_discrete/6_aorta_remeshed_pt5mm_2cm_extender_layers_constriction.vtu . 

# extracts relevant portion of mesh 
python3 ~/copies_scripts/remove_unnecessary_eulerian_space.py $TOTAL_TASKS

# run integral metrics 
sbatch ~/copies_scripts/run_integral_metrics.sh

# adds face data 
python3 ~/copies_scripts/add_faces.py

python3 ~/copies_scripts/fix_pvd_files.py

python3 ~/copies_scripts/convert_vtu_vertices_to_csv.py

module load matlab/R2017b
cp ../aortic_no_partition*final_data.mat . 
matlab -nodesktop -nodisplay -r 'addpath ~/valve_generator; generate_cells_file; "matlab generate cells complete"; exit;'

python3 ~/copies_scripts/convert_csv_with_cells.py


SESSION_NAME_PARAVIEW="~/copies_scripts/bicuspid_slices_paraview.py"
# SESSION_NAME_PARAVIEW="~/copies_scripts/bicuspid_slices_paraview_paper.py"

# renders movie 
# assumes located in main sim directory 
cd .. 
python ~/copies_scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS

SESSION_NAME_PARAVIEW_PAPER="~/copies_scripts/bicuspid_slices_paraview_paper.py"
python ~/copies_scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_PAPER $TOTAL_TASKS

SESSION_NAME_PARAVIEW_TOP="~/copies_scripts/top_view_valve_0_paper.py"
python ~/copies_scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_TOP $TOTAL_TASKS

SESSION_NAME_PARAVIEW_VERTICAL="~/copies_scripts/bicuspid_slices_paraview_vertical.py"
python ~/copies_scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_VERTICAL $TOTAL_TASKS