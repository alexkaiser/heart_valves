#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=hist3_post
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --partition=amarsden


module purge
module load gcc/8.1.0
module load openmpi/2.0.2

source ~/.bash_profile

module load python/3.9.0
module load matlab/R2017b


BOUNDARY_MESH_NAME=1_aorta_lv_ross_1_extenders_remesh_pt25mm_cap_CM_UNITS.stl
VESSEL_BASE_NAME=aorta_ross_1_pt5mm_192
VESSEL_WITH_FACES_NAME=aorta_ross_1_pt5mm_192.vtu

if [ -e $BOUNDARY_MESH_NAME ]
then
    echo $BOUNDARY_MESH_NAME " found"
else
    cp ~/ross1_meshes/$BOUNDARY_MESH_NAME . 
fi

if [ -e $VESSEL_WITH_FACES_NAME ]
then
    echo $VESSEL_WITH_FACES_NAME " found"
else
    cp ~/ross1_meshes/$VESSEL_WITH_FACES_NAME . 
fi


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


# extracts relevant portion of mesh 
python3 ~/heart_valves/scripts/remove_unnecessary_eulerian_space.py $TOTAL_TASKS $BOUNDARY_MESH_NAME

wait

# run integral metrics 
# sbatch ~/heart_valves/scripts/run_integral_metrics.sh
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py annulus_normal_projected & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_3 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_6 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_9 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_12 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/heart_valves/scripts/integrals_contour_12_2.py contour_21 & 


# adds face data 
python3 ~/heart_valves/scripts/add_faces.py $VESSEL_BASE_NAME $VESSEL_WITH_FACES_NAME

python3 ~/heart_valves/scripts/fix_pvd_files.py

python3 ~/heart_valves/scripts/convert_vtu_vertices_to_csv.py

cp ../bc_data.mat . 
cp ../aortic_no_partition*final_data.mat . 
matlab -nodesktop -nodisplay -r 'addpath ~/heart_valves/valve_generator; generate_cells_file; "matlab generate cells complete"; exit;'

python3 ~/heart_valves/scripts/convert_csv_with_cells.py

# ensure integral metrics have finished 
wait 

# renders movie 
# assumes located in main sim directory 
cd .. 

PATH=/home/groups/amarsden/ParaView-5.13.3-osmesa-MPI-Linux-Python3.10-x86_64/bin:$PATH

SESSION_NAME_PARAVIEW="~/heart_valves/scripts/ross_1_multiview.py"
python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW $TOTAL_TASKS

# SESSION_NAME_PARAVIEW_TOP="~/heart_valves/scripts/top_view_valve_0_paper.py"
# python ~/heart_valves/scripts/run_parallel_movie_paraview.py $SESSION_NAME_PARAVIEW_TOP $TOTAL_TASKS


