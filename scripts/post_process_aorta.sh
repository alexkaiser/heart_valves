#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=8:00:00
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

cp ~/mitral_fully_discrete/2_aorta_remeshed_pt5mm_capped.vtp . 
cp ~/mitral_fully_discrete/4_aorta_remeshed_pt25mm_3cm_extender_layers_constriction.vtu . 
cp ~/mitral_fully_discrete/6_aorta_remeshed_pt5mm_2cm_extender_layers_constriction.vtu . 

# extracts relevant portion of mesh 
python3 ~/copies_scripts/remove_unnecessary_eulerian_space.py $SLURM_NTASKS

# adds face data 
python3 ~/copies_scripts/add_faces.py

SESSION_NAME_PARAVIEW="~/copies_scripts/bicuspid_slices_paraview.py"
# SESSION_NAME_PARAVIEW="~/copies_scripts/bicuspid_slices_paraview_paper.py"

# renders movie 
# assumes located in main sim directory 
cd .. 
python ~/copies_scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $SLURM_NTASKS
