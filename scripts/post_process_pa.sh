#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mem=185GB
#SBATCH --job-name=post_process
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
# #SBATCH --partition=willhies
#SBATCH --partition=willhies,amarsden



# load stuff for movie making 
source ~/.bash_profile

# change to viz directory 
cd viz_IB3d_tree_cycle_256 

# average in time 
python3 ~/heart_valves/scripts/run_temporal_averaging.py
wait 

# change to the averaged directory 
cd vis_data_averaged_cycle_2

# pull data to interpolate 
cp ~/HealthyNative_vtk/*vtk . 

# run interpolation 
python3 ~/heart_valves/scripts/interpolate_mesh.py

# add faces to vessel scripts 
python3 ~/heart_valves/scripts/add_faces.py

# finally render 
pvbatch ~/heart_valves/scripts/axial_view_points.py 0 &
pvbatch ~/heart_valves/scripts/axial_view_points.py 1 &
pvbatch ~/heart_valves/scripts/axial_view_points.py 2 &
pvbatch ~/heart_valves/scripts/differences_multiview_points.py &
pvbatch ~/heart_valves/scripts/sagittal_fine_cropped_points.py &

# and integral metrics 
# pvbatch ~/heart_valves/scripts/convergence_study_integrals.py &
# pvbatch ~/heart_valves/scripts/integral_metrics.py 0 &
# pvbatch ~/heart_valves/scripts/integral_metrics.py 1 &
# pvbatch ~/heart_valves/scripts/integral_metrics.py 2 &

wait 
