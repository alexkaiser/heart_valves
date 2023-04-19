#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --time=23:59:00
#SBATCH --mem=60GB
#SBATCH --job-name=integral_metrics
#SBATCH --mail-user=adkaiser@gmail.com
#SBATCH --mail-type=ALL
# #SBATCH --partition=amarsden
#SBATCH --partition=willhies,amarsden
# #SBATCH --exclude=sh-107-[59-64]

module purge
module load gcc/8.1.0
module load openmpi/2.0.2

source ~/.bash_profile

cd viz_IB3d_tree_cycle_256

/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py annulus_normal_projected & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_3 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_6 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_9 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_12 & 
/home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_21 & 


# cd aortic_9952752_384_7da89c7_final_setup_r_5pt6_p0_94_plus_15_sys/viz_IB3d_tree_cycle_256

# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py annulus_normal_projected & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_3 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_6 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_9 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_12 & 

# cd ../..

# cd aortic_9998511_384_ef4cbc_final_setup_2_aa5e813_mesh_comm_0_fused/viz_IB3d_tree_cycle_256

# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py annulus_normal_projected & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_3 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_6 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_9 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_12 & 

# cd ../..

# cd aortic_9999320_384_ef4cbc_final_setup_2_d6108c3_mesh_comm_1_fused/viz_IB3d_tree_cycle_256

# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py annulus_normal_projected & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_3 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_6 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_9 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_12 & 

# cd ../..

# cd aortic_9999816_384_ef4cbc_final_setup_2_ee251fa_mesh_comm_2_fused/viz_IB3d_tree_cycle_256

# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py annulus_normal_projected & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_3 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_6 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_9 & 
# /home/groups/amarsden/ParaView-5.9.0-osmesa-MPI-Linux-Python3.8-64bit/bin/pvbatch ~/copies_scripts/integrals_contour_12_2.py contour_12 & 

# cd ../..


wait 
echo "integral_metrics done"



