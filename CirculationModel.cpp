// Filename: CirculationModel.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified by Alex Kaiser, 7/2016

#include "CirculationModel.h"

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
#include <cmath>

#include <Eigen/Dense>
using namespace Eigen;

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

#define MIN_PER_L_T0_SEC_PER_ML 60.0e-3   // 60/1000
#define MMHG_TO_CGS 1333.22368


#define USE_WINDKESSEL
// If not defined then this only computes fluxes
// Pressure is set to zero
// Currently hard coded for upper boundary

namespace
{
    // Name of output file.
    static const string DATA_FILE_NAME = "bc_data.m";

    // constants 
    static const double C_PA =  4.12; // Pulmonary artery compliance, ml / mmHg 
    static const double C_PV = 10.0;  // Pulmonary vein compliance, ml / mmHg 
    static const double C_LA =  1.6;  // Left atrial compliance ml / mmHg  
    
    static const double R_P  = (9.0/5.6) * MIN_PER_L_T0_SEC_PER_ML; // Pulmonary resistance, mmHg / (ml/s)
    
    static const double beat_time = 0.8; 
    static const double T_on = .53;   // Pulmonary valve open 
    static const double T_off = .75;  // Pulmonary valve closes  
    static const double T_peak = T_on + 0.4 * (T_off - T_on); // Peak pulmonary valve flow 
    static const double stroke_volume = 75; // ml 
    static const double h = 2.0 * stroke_volume / (T_off - T_on); // Peak flow to get given stroke volume
    
    
    inline double compute_Q_R(double t){
        // Triangle wave flux 
        
        double t_reduced = t - beat_time * floor(t/beat_time); 
        
        if (t_reduced <= T_on)
            return 0.0; 
        else if (t_reduced <= T_peak)
            return ( h/(T_peak - T_on) )*t_reduced - (h/(T_peak - T_on)  )*T_on;
        else if (t_reduced <= T_off)
            return (-h/(T_off - T_peak))*t_reduced + (h/(T_off -  T_peak))*T_off;
        else if (t_reduced <= beat_time)
            return 0.0;
        
        TBOX_ERROR("Valid time for flux not found.");
        return 0.0;
    }
    
    // Backward Euler update for windkessel model.
    inline void
    windkessel_be_update(double& Q_R, double& P_PA, double& Q_P, double& P_LA, const double Q_mi, const double t, const double dt)
    {
     
        Q_R = compute_Q_R(t);
        
        double a = C_PA/dt + 1/R_P; 
        double b = -1/R_P; 
        double c = -1/R_P; 
        double d = (C_PV + C_LA)/dt + 1/R_P; 
        
        double rhs[2];  
        rhs[0] = (C_PA/dt)*P_PA + Q_R;
        rhs[1] = ((C_PV + C_LA)/dt)*P_LA - Q_mi;
        
        double det = a*d - b*c; 
        
        // Closed form linear system solution 
        P_PA = (1/det) * ( d*rhs[0] + -b*rhs[1]);
        P_LA = (1/det) * (-c*rhs[0] +  a*rhs[1]);
                
        Q_P = (1/R_P) * (P_PA - P_LA);
                
        return;
    } // windkessel_be_update

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel::CirculationModel(const string& object_name, double P_PA_0, double P_LA_0, double t, bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_time(t),
      d_nsrc(1),           // number of sets of variables
      d_psrc(d_nsrc, 0.0), // pressure
      d_qsrc(d_nsrc, 0.0), // flux
      d_srcname(d_nsrc),
      d_P_PA(P_PA_0),
      d_P_LA(P_LA_0), 
      d_Q_R(0.0),
      d_Q_P(0.0), 
      d_Q_mi(0.0),
      d_bdry_interface_level_number(numeric_limits<int>::max())
{
#if !defined(NDEBUG)
    assert(!object_name.empty());
#endif
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
    else
    {
        //   nsrcs = the number of sources in the valve tester:
        //           (1) left atrium
        d_srcname[0] = "left atrium       ";        
    }
    return;
} // CirculationModel

CirculationModel::~CirculationModel()
{
    return;
} // ~CirculationModel

void
CirculationModel::advanceTimeDependentData(const double dt,
                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const int U_idx,
                                           const int /*P_idx*/,
                                           const int /*wgt_cc_idx*/,
                                           const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    std::fill(d_qsrc.begin(), d_qsrc.end(), 0.0);
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
                // const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }

                static const int axis = 2;  // Always z axis here
                const int side = 1;         // Compute flux at the top only

                const bool is_lower = side == 0;
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
                        
                        // no conditional here, just add the flux in
                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                        if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                        {
                            double dA = n[axis] * dV / dx[axis];
                            d_qsrc[0] += (*U_data)(i_s)*dA;
                            
                            // pout << "adding " << (*U_data)(i_s)*dA << "to the flux\n";  
                            
                        }
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&d_qsrc[0], d_nsrc);

    // pout << "computed flux = " << d_qsrc[0] << "\n"; 
    
    
    #ifdef USE_WINDKESSEL
        // The downstream (Atrial) pressure is determined by a zero-d model 
        const double t = d_time;
    
        double& Q_R  = d_Q_R;
        double& P_PA = d_P_PA;
        double& Q_P  = d_Q_P;
        double& P_LA = d_P_LA;
    
        // Mitral flux is inward flow, d_qsrc is outward flux
        d_Q_mi = -d_qsrc[0];
        
        windkessel_be_update(Q_R, P_PA, Q_P, P_LA, d_Q_mi, t, dt);
        
        // model in mmHg, to CGS for solver pressure
        d_psrc[0] = d_P_LA * MMHG_TO_CGS;
    #else 
        // Downstream pressure is fixed to zero 
        d_psrc[0] = 0.0;
    #endif

    // Update the current time.
    d_time += dt;

    #ifdef USE_WINDKESSEL

        // Output the updated values.
        const long precision = plog.precision();
        plog.unsetf(ios_base::showpos);
        plog.unsetf(ios_base::scientific);

        plog.precision(12);

        plog << "============================================================================\n"
             << "Circulation model variables at time " << d_time << ":\n";

        plog << "P_PA (mmHg)\t P_LA (mmHg)\t Q_R (ml/s)\t Q_P (ml/s)\t Q_mi (ml/s)\n";
        plog.setf(ios_base::showpos);
        plog.setf(ios_base::scientific);
        
        plog << d_P_PA << ",\t " << d_P_LA << ",\t " << d_Q_R << ",\t " << d_Q_P << ",\t " << d_Q_mi << "\n";
        plog << "============================================================================\n";

        plog.unsetf(ios_base::showpos);
        plog.unsetf(ios_base::scientific);
        plog.precision(precision);
    
    #endif
    // Write the current state to disk.
    writeDataFile();
    return;
} // advanceTimeDependentData

void
CirculationModel::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_time", d_time);
    db->putInteger("d_nsrc", d_nsrc);
    db->putDoubleArray("d_qsrc", &d_qsrc[0], d_nsrc);
    db->putDoubleArray("d_psrc", &d_psrc[0], d_nsrc);
    db->putStringArray("d_srcname", &d_srcname[0], d_nsrc);
    db->putDouble("d_P_PA", d_P_PA);
    db->putDouble("d_P_LA", d_P_LA);
    db->putDouble("d_Q_R", d_Q_R);
    db->putDouble("d_Q_P", d_Q_P);
    db->putDouble("d_Q_mi", d_Q_mi);
    db->putInteger("d_bdry_interface_level_number", d_bdry_interface_level_number);
    return;
} // putToDatabase



/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CirculationModel::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time \t P_PA (mmHg)\t d_P_LA (mmHg)\t Q_R (ml/s)\t Q_P (ml/s)\t Q_mi (ml/s)"
                 << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        for (int n = 0; n < d_nsrc; ++n)
        {
            fout << d_time;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(10);
            fout << " " << d_P_PA << " " << d_P_LA << " " << d_Q_R << " " << d_Q_P << " " << d_Q_mi << "; \n";
        }
    }
    return;
} // writeDataFile

void
CirculationModel::getFromRestart()
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

    d_time = db->getDouble("d_time");
    d_nsrc = db->getInteger("d_nsrc");
    d_qsrc.resize(d_nsrc);
    d_psrc.resize(d_nsrc);
    d_srcname.resize(d_nsrc);
    db->getDoubleArray("d_qsrc", &d_qsrc[0], d_nsrc);
    db->getDoubleArray("d_psrc", &d_psrc[0], d_nsrc);
    db->getStringArray("d_srcname", &d_srcname[0], d_nsrc);
    d_P_PA = db->getDouble("d_P_PA");
    d_P_LA = db->getDouble("d_P_LA");
    d_Q_R  = db->getDouble("d_Q_R");
    d_Q_P  = db->getDouble("d_Q_P");
    d_Q_mi = db->getDouble("d_Q_mi");
    d_bdry_interface_level_number = db->getInteger("d_bdry_interface_level_number");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////