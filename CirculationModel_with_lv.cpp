// Filename: CirculationModel_with_lv.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_with_lv.h"

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

CirculationModel_with_lv::CirculationModel_with_lv(const fourier_series_data *fourier_aorta, 
                                                   const fourier_series_data *fourier_atrium, 
                                                   const fourier_series_data *fourier_ventricle, 
                                                   const double  radius_aorta,
                                                   const double  radius_atrium,
                                                   const double* center_aorta,
                                                   const double* center_atrium, 
                                                   const double  cycle_duration,
                                                   const double  t_offset_bcs_unscaled, 
                                                   const double  initial_time)
    : d_fourier_aorta (fourier_aorta), 
      d_fourier_atrium(fourier_atrium),       
      d_fourier_ventricle(fourier_ventricle), 
      d_radius_aorta  (radius_aorta),
      d_radius_atrium (radius_atrium),
      d_center_aorta  (center_aorta),
      d_center_atrium (center_atrium),       
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_current_idx_series(0),
      d_Q_aorta      (0.0), 
      d_Q_left_atrium(0.0),
      d_Q_mitral     (0.0),
      d_time(initial_time),
      d_object_name("circ_model_with_lv")  // constant name here 
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
    
    return;
} // CirculationModel

CirculationModel_with_lv::~CirculationModel_with_lv()
{
    return;
} // ~CirculationModel_with_lv


void CirculationModel_with_lv::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    
    double Q_aorta_local = 0.0; 
    double Q_left_atrium_local = 0.0; 

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
                double X[NDIM];

                // always looking for z flux here 
                // side is always 1, top of box 
                static const int axis = 2;
                int side = 1; 
                const bool is_lower = (side == 0);
                if (pgeom->getTouchesRegularBoundary(axis, side))
                {

//                    radius and position of sources are handled by input data 
//                    const double rsrc = d_rsrc[side];
//                    const Point& posn = d_posn[side];
                    
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
                        double dist_sq_aorta = 0.0;
                        double dist_sq_atrium = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                            if (d != axis){
                                dist_sq_aorta  += pow(X[d] - d_center_aorta[d],  2.0);
                                dist_sq_atrium += pow(X[d] - d_center_atrium[d], 2.0);
                            }
                        }
                        const double dist_aorta  = sqrt(dist_sq_aorta);
                        const double dist_atrium = sqrt(dist_sq_atrium);

                        if (dist_aorta < d_radius_aorta)
                        {
                            const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                            if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                            {
                                double dA = n[axis] * dV / dx[axis];
                                Q_aorta_local += (*U_data)(i_s)*dA;
                            }
                        }

                        if (dist_atrium < d_radius_atrium)
                        {
                            const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                            if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                            {
                                double dA = n[axis] * dV / dx[axis];
                                Q_left_atrium_local += (*U_data)(i_s)*dA;
                            }
                        }
                    }
                }
                
            }
        }
    }

    d_Q_aorta       = SAMRAI_MPI::sumReduction(Q_aorta_local);
    d_Q_left_atrium = SAMRAI_MPI::sumReduction(Q_left_atrium_local);

    d_time += dt; 

    // compute which index in the Fourier series we need here 
    // always use a time in current cycle 
    double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration); 

    // fourier series has its own period, scale to that 
    double t_scaled = t_reduced * (d_fourier_aorta->L  / d_cycle_duration); 

    // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
    double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

    // Fourier data here
    // index without periodicity 
    unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_aorta->dt));
    
    // // take periodic reduction
    unsigned int d_current_idx_series = k % (d_fourier_aorta->N_times);

    writeDataFile(); 

} // advanceTimeDependentData

void CirculationModel_with_lv::set_Q_mitral(double Q_mitral){
    d_Q_mitral = Q_mitral; 
}


void
CirculationModel_with_lv::putToDatabase(Pointer<Database> db)
{

    db->putInteger("d_current_idx_series", d_current_idx_series); 
    db->putDouble("d_Q_aorta", d_Q_aorta); 
    db->putDouble("d_Q_left_atrium", d_Q_left_atrium);
    db->putDouble("d_Q_mitral", d_Q_mitral);
    db->putDouble("d_time", d_time); 
    return; 
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CirculationModel_with_lv::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time \t P_aorta (mmHg)\t P_atrium (mmHg)\t P_ventricle (mmHg)\t Q_P (ml/s)\t Q_mi (ml/s)"
                 << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);

        fout << d_time;
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        double P_aorta     = MMHG_TO_CGS * d_fourier_aorta->values[d_current_idx_series]; 
        double P_atrium    = MMHG_TO_CGS * d_fourier_atrium->values[d_current_idx_series];
        double P_ventricle = MMHG_TO_CGS * d_fourier_ventricle->values[d_current_idx_series];
        fout << " " << P_aorta <<  " " << P_atrium << " " << P_ventricle << " " << d_Q_aorta << " " << d_Q_left_atrium << " " << d_Q_mitral << "; \n";

    }

    return;
} // writeDataFile

void
CirculationModel_with_lv::getFromRestart()
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
    d_Q_aorta            = db->getDouble("d_Q_aorta"); 
    d_Q_left_atrium      = db->getDouble("d_Q_left_atrium");
    d_Q_mitral           = db->getDouble("d_Q_mitral");
    d_time               = db->getDouble("d_time");

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////