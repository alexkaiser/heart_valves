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
module load matlab/R2017b

# generate bc data 
if [ -e bc_data.mat ]
then
    echo "bc_data.mat found"
else
    matlab -nodesktop -nodisplay -r 'addpath ~/heart_valves/valve_generator; bc_data; exit;' & 
fi


# change dir even if not needed 
cd viz_IB3d_tree_cycle_256 

TOTAL_TASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))

cp ~/heart_valves/2_aorta_remeshed_pt5mm_capped.vtp . 

# extracts relevant portion of mesh 
python3 ~/heart_valves/scripts/remove_unnecessary_eulerian_space.py $TOTAL_TASKS


wait

# run integral metrics 
# sbatch ~/heart_valves/scripts/run_integral_metrics.sh
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py annulus_normal_projected & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_3 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_6 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_9 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_12 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_21 & 


# adds face data 
python3 ~/heart_valves/scripts/add_faces.py

python3 ~/heart_valves/scripts/fix_pvd_files.py

python3 ~/heart_valves/scripts/convert_vtu_vertices_to_csv.py

cp ../bc_data.mat . 
cp ../aortic_no_partition*final_data.mat . 
matlab -nodesktop -nodisplay -r 'addpath ~/heart_valves/valve_generator; generate_cells_file; "matlab generate cells complete"; exit;'

python3 ~/heart_valves/scripts/convert_csv_with_cells.py

# ensure integral metrics have finished 
wait 

matlab -nodesktop -nodisplay -r 'addpath ~/heart_valves/valve_generator; run_shape_analysis_local; exit;'

# SESSION_NAME_PARAVIEW="~/heart_valves/scripts/bicuspid_slices_paraview_paper.py"

# renders movie 
# assumes located in main sim directory 
cd .. 

SESSION_NAME_PARAVIEW="~/heart_valves/scripts/bicuspid_slices_paraview.py"
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS

# SESSION_NAME_PARAVIEW_PAPER="~/heart_valves/scripts/bicuspid_slices_paraview_paper.py"
# python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_PAPER $TOTAL_TASKS

SESSION_NAME_PARAVIEW_TOP="~/heart_valves/scripts/top_view_valve_0_paper.py"
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_TOP $TOTAL_TASKS

SESSION_NAME_PARAVIEW_VERTICAL="~/heart_valves/scripts/bicuspid_slices_paraview_vertical.py"
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_VERTICAL $TOTAL_TASKS
