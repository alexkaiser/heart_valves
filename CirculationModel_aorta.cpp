// Filename: CirculationModel_aorta.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_aorta.h"
#include "pnpoly.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

#include <Eigen/Dense>
using namespace Eigen;

namespace
{
// Name of output file.
static const string DATA_FILE_NAME = "bc_data.m";

} 

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel_aorta::CirculationModel_aorta(Pointer<Database> input_db, 
                                               const fourier_series_data *fourier_ventricle, 
                                               string ventricle_vertices_file_name,
                                               string aorta_vertices_file_name,
                                               const double  cycle_duration,
                                               const double  t_offset_bcs_unscaled, 
                                               const double  initial_time, 
                                               double P_initial_aorta,
                                               bool rcr_bcs_on)
    : 
      d_object_name("circ_model_aorta"),  // constant name here  
      d_registered_for_restart(true),      // always true
      d_fourier_ventricle(fourier_ventricle), 
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_current_idx_series(0),
      d_Q_ventricle(0.0), 
      d_Q_aorta(0.0),
      d_time(initial_time), 
      d_aorta_P(P_initial_aorta), 
      d_aorta_P_Wk(P_initial_aorta),
      d_area_ventricle(0.0),
      d_area_aorta(0.0),
      d_area_initialized(false), 
      d_rcr_bcs_on(rcr_bcs_on)
{
    
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    
    if (d_rcr_bcs_on){
        if (input_db){
            // left and right equal for now
            d_aorta_R_proximal = input_db->getDouble("R_proximal");
            d_aorta_R_distal   = input_db->getDouble("R_distal");
            d_aorta_C          = input_db->getDouble("C");

            std::cout << "input db got values:\n";
            std::cout << "input db got values R_proximal = " << d_aorta_R_proximal << "\tR_distal = " << d_aorta_R_distal << "\tC = " << d_aorta_C << "\n";   
        }
        else {
            TBOX_ERROR("Must provide valid input_db");
        }
    }

    double x,x_prev,y,y_prev,z,z_prev; 
    double tol = 1.0e-2; 

    // read vertices from file 
    ifstream ventricle_file(ventricle_vertices_file_name.c_str(), ios::in);

    if(!ventricle_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    ventricle_file >> d_n_pts_ventricle; 
    
    d_ventricle_points_idx1 = new double[d_n_pts_ventricle]; 
    d_ventricle_points_idx2 = new double[d_n_pts_ventricle]; 

    for (int i=0; i<d_n_pts_ventricle; i++){
        ventricle_file >> x; 
        ventricle_file >> d_ventricle_points_idx1[i]; 
        ventricle_file >> d_ventricle_points_idx2[i];
        
        if (i>0){
            if (fabs(x_prev - x) > tol){
                TBOX_ERROR("x coordinates must be consistent\n"); 
            }
        }
        x_prev = x; 

    }
    pout << "to ventricle file close\n"; 
    ventricle_file.close(); 
    // hardcode to top 
    d_ventricle_axis = 0; 
    d_ventricle_side = 1; 

    // read vertices from file 
    ifstream aorta_file(aorta_vertices_file_name.c_str(), ios::in);

    if(!aorta_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    aorta_file >> d_n_pts_aorta; 
    
    d_aorta_points_idx1 = new double[d_n_pts_aorta]; 
    d_aorta_points_idx2 = new double[d_n_pts_aorta]; 

    for (int i=0; i<d_n_pts_aorta; i++){
        aorta_file >> d_aorta_points_idx1[i]; 
        aorta_file >> d_aorta_points_idx2[i]; 
        aorta_file >> z;
        

        if (i>0){
            if (fabs(z_prev - z) > tol){
                TBOX_ERROR("z coordinates must be consistent\n"); 
            }
        }
        z_prev = z; 

    }
    pout << "to aorta file close\n"; 
    aorta_file.close(); 
    d_aorta_axis = 2; 
    d_aorta_side = 1; 

    pout << "passed contstructor\n"; 

    pout << "initial aorta pressure = " << P_initial_aorta << ", P_wk = " << d_aorta_P << "\n"; 

    return;
} // CirculationModel

CirculationModel_aorta::~CirculationModel_aorta()
{
    return;
} // ~CirculationModel_aorta


void CirculationModel_aorta::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    
    double Q_ventricle_local = 0.0; 
    double Q_aorta_local = 0.0; 

    double area_ventricle_local = 0.0; 
    double area_aorta_local = 0.0; 

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }

                for(int axis=0; axis<3; axis++)
                {
                    for(int side=0; side<2; side++)
                    {
                        const bool is_lower = (side == 0);
                        if (pgeom->getTouchesRegularBoundary(axis, side))
                        {
                            
                            Vector n;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                            }
                            Box<NDIM> side_box = patch_box;
                            if (is_lower)
                            {
                                side_box.lower(axis) = patch_box.lower(axis);
                                side_box.upper(axis) = patch_box.lower(axis);
                            }
                            else
                            {
                                side_box.lower(axis) = patch_box.upper(axis) + 1;
                                side_box.upper(axis) = patch_box.upper(axis) + 1;
                            }
                            for (Box<NDIM>::Iterator b(side_box); b; b++)
                            {
                                const Index<NDIM>& i = b();

                                double X[NDIM];
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                }

                                double X_in_plane_1 = 0.0; 
                                double X_in_plane_2 = 0.0; 
                                if (axis == 0)
                                {
                                    X_in_plane_1 = X[1]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 1)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 2)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[1]; 
                                }
                                else{
                                    TBOX_ERROR("Invalid value of axis\n"); 
                                }

                                const int in_ventricle  = this->point_in_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_aorta      = this->point_in_aorta    (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_ventricle && in_aorta){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }

                                if (in_ventricle)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_ventricle_local += (*U_data)(i_s)* n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_ventricle_local += dA;
                                        }

                                    }
                                }

                                if (in_aorta)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_aorta_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_aorta_local += dA;
                                        }

                                    }
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    d_Q_ventricle = SAMRAI_MPI::sumReduction(Q_ventricle_local);
    d_Q_aorta        = SAMRAI_MPI::sumReduction(Q_aorta_local);

    if (!d_area_initialized){
        d_area_ventricle = SAMRAI_MPI::sumReduction(area_ventricle_local);
        d_area_aorta        = SAMRAI_MPI::sumReduction(area_aorta_local);  
        d_area_initialized = true;       
    }

    if (d_rcr_bcs_on){
        // The downstream pressure is determined by a three-element Windkessel model.

        d_aorta_P_Wk = ((d_aorta_C / dt) * d_aorta_P_Wk + d_Q_aorta) / (d_aorta_C / dt + 1.0 / d_aorta_R_distal);        
        d_aorta_P = d_aorta_P_Wk + d_aorta_R_proximal * d_Q_aorta;

    }

    // print_summary();

    // bool debug_out_areas = false; 
    // if (debug_out_areas){
    //     pout << "d_area_ventricle = " << d_area_ventricle << "\n"; 
    //     pout << "d_area_aorta = " << d_area_aorta << "\n"; 
    // }

    d_time += dt; 

    // compute which index in the Fourier series we need here 
    // always use a time in current cycle 
    double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration); 

    // fourier series has its own period, scale to that 
    double t_scaled = t_reduced * (d_fourier_ventricle->L  / d_cycle_duration); 

    // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
    double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

    // Fourier data here
    // index without periodicity 
    unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_ventricle->dt));
    
    // // take periodic reduction
    d_current_idx_series = k % (d_fourier_ventricle->N_times);

    // bool debug_out = false; 
    // if (debug_out){
    //     pout << "circ mode: d_time = " << d_time << ", d_current_idx_series = " << d_current_idx_series << "\n"; 
    //     pout << "t_reduced = " << t_reduced << " t_scaled = " << t_scaled << " t_scaled_offset = " << t_scaled_offset << "\n"; 
    //     pout << "k (unreduced idx) = " << k << " d_current_idx_series = " << d_current_idx_series << "\n\n"; 
    // }


    writeDataFile(); 

} // advanceTimeDependentData

void CirculationModel_aorta::set_Q_valve(double Q_valve){
    d_Q_valve = Q_valve; 
}



void
CirculationModel_aorta::putToDatabase(Pointer<Database> db)
{

    db->putInteger("d_current_idx_series", d_current_idx_series); 
    db->putDouble("d_Q_ventricle", d_Q_ventricle); 
    db->putDouble("d_Q_aorta", d_Q_aorta);
    db->putDouble("d_Q_valve", d_Q_valve);
    db->putDouble("d_aorta_P", d_aorta_P);
    db->putDouble("d_aorta_P_Wk", d_aorta_P_Wk);
    db->putDouble("d_time", d_time); 
    db->putBool("d_rcr_bcs_on", d_rcr_bcs_on); 
    return; 
} // putToDatabase

void CirculationModel_aorta::print_summary(){

    double P_ventricle = d_fourier_ventricle->values[d_current_idx_series]; 
    double P_aorta; 

    if (d_rcr_bcs_on){
        P_aorta        = d_aorta_P / MMHG_TO_CGS;
    }
    else{
        TBOX_ERROR("Not implemented\n"); 
        // P_aorta        = d_fourier_aorta->values[d_current_idx_series];
    }

    pout << "rcr_bcs_on = " << d_rcr_bcs_on << "\n"; 
    pout << "% time \t       P_ventricle (mmHg)\t   P_aorta (mmHg)\t  Q_ventricle (ml/s)\t    d_Q_aorta (ml/s)\t  d_Q_valve (ml/s)\t  Q_current_idx_series \t idx" ;
    if (d_rcr_bcs_on){
        pout << "\t aorta_P_Wk "; 
    }
    pout << "\n";
    pout << d_time << " " << P_ventricle <<  " " << P_aorta << " " << d_Q_ventricle << " " << d_Q_aorta << " " << d_Q_valve << " " << d_current_idx_series; 
    if (d_rcr_bcs_on){
        pout  << " " << d_aorta_P_Wk; 
    }
    pout << "\n";

}

int CirculationModel_aorta::point_in_ventricle(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_ventricle_axis) || (side != d_ventricle_side))
        return 0; 

    return pnpoly(d_n_pts_ventricle, d_ventricle_points_idx1, d_ventricle_points_idx2, testx, testy); 
}

int CirculationModel_aorta::point_in_aorta(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_aorta_axis) || (side != d_aorta_side))
        return 0; 

    return pnpoly(d_n_pts_aorta, d_aorta_points_idx1, d_aorta_points_idx2, testx, testy); 
}


void CirculationModel_aorta::write_plot_code()
{

    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);
        fout << "];\n";  
        fout << "times       =  bc_vals(:,1);\n"; 
        fout << "p_lv        =  bc_vals(:,2);\n"; 
        fout << "p_aorta     =  bc_vals(:,3); \n"; 
        fout << "q_ventricle = -bc_vals(:,4);\n"; 
        fout << "q_aorta     =  bc_vals(:,5);\n"; 
        fout << "q_valve     =  bc_vals(:,6);\n"; 
        fout << "p_wk        =  bc_vals(:,7);\n"; 
        fout << "subplot(2,1,1)\n"; 
        fout << "plot(times, p_aorta, 'k')\n"; 
        fout << "hold on\n"; 
        fout << "plot(times, p_wk, ':k')\n"; 
        fout << "plot(times, p_lv, '--k')\n"; 
        fout << "legend('P_{Ao}', 'P_{Wk}', 'P_{LV}', 'Location','NorthEastOutside');\n"; 
        fout << "xlabel('t (s)');\n"; 
        fout << "ylabel('P (mmHg)');\n"; 
        fout << "subplot(2,1,2)\n"; 
        fout << "plot(times, q_aorta, 'k')\n"; 
        fout << "hold on\n"; 
        fout << "dt = times(2,1) - times(1);\n"; 
        fout << "net_flux = dt*cumsum(q_aorta);\n"; 
        fout << "plot(times, net_flux, '--k')\n"; 
        fout << "plot(times, q_ventricle)\n"; 
        fout << "plot(times, q_valve)\n"; 
        fout << "plot(bc_vals(:,1), 0*net_flux, ':k')\n"; 
        fout << "legend('Q', 'net Q', 'Q ventricle', 'Q valve', 'Location','NorthEastOutside')\n"; 
        fout << "xlabel('t (s)')\n"; 
        fout << "ylabel('Flow (ml/s), Net Flow (ml)')\n"; 
        fout << "set(fig, 'Position', [100, 100, 1000, 750])\n"; 
        fout << "set(fig,'PaperPositionMode','auto')\n"; 
        fout << "printfig(fig, 'bc_model_variables')\n"; 
        fout << "min_p_aorta_after_first_beat = min(p_aorta(floor(end/3):end))\n"; 
        fout << "max_p_aorta_after_first_beat = max(p_aorta(floor(end/3):end))\n"; 
        fout << "mean_p_aorta = mean(p_aorta)\n"; 
        fout << "mean_p_wk    = mean(p_wk)\n"; 
        fout << "mean_p_lv    = mean(p_lv)\n"; 

    }
    return;

}



/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
    CirculationModel_aorta::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time \t P_ventricle (mmHg)\t P_aorta (mmHg)\t d_Q_ventricle (ml/s)\t d_Q_aorta (ml/s)\t d_Q_valve (ml/s)"
                 << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);

        fout << d_time;
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        double P_ventricle = d_fourier_ventricle->values[d_current_idx_series]; 

        double P_aorta = 0.0; 

        if (d_rcr_bcs_on){
            P_aorta        = d_aorta_P/MMHG_TO_CGS;
        }
        else{
            TBOX_ERROR("not implemented\n"); 
            // P_aorta        = d_fourier_aorta->values[d_current_idx_series];
        }

        fout << " " << P_ventricle <<  " " << P_aorta;
        fout << " " << d_Q_ventricle << " " << d_Q_aorta << " " << d_Q_valve;         
        fout << " " << d_aorta_P_Wk/MMHG_TO_CGS;
        fout << "; \n";

    }

    return;
} // writeDataFile

void
CirculationModel_aorta::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_current_idx_series = db->getInteger("d_current_idx_series"); 
    d_Q_ventricle        = db->getDouble("d_Q_ventricle"); 
    d_Q_aorta            = db->getDouble("d_Q_aorta");
    d_Q_valve            = db->getDouble("d_Q_valve");
    d_aorta_P            = db->getDouble("d_aorta_P");
    d_aorta_P_Wk         = db->getDouble("d_aorta_P_Wk");
    d_time               = db->getDouble("d_time");
    d_rcr_bcs_on         = db->getBool("d_rcr_bcs_on"); 
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////