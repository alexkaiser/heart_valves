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
    static const string DATA_FILE_NAME = "bc_data.txt";

    // Three-element windkessel model
    // Both resistances set to 0.01 mmHg s / L
    // Taken from Peskin's book as fwd valve resistances which do not do much
    // Compliance is pulmonary artery complaince
    static double R_P = 0.01 * MIN_PER_L_T0_SEC_PER_ML; // peripheral resistance (mmHg ml^-1 s)
    static double R_C = 0.01 * MIN_PER_L_T0_SEC_PER_ML; // characteristic resistance (mmHg ml^-1 s)
    static double C   = 0.01 * 1e3;   // pulmonary vein compliance, ml / mmHg

    // Time required to "ramp up" the pressure in the Windkessel model.
    static double P_ramp = 5.0;
    static double t_ramp = 0.05;

    // Backward Euler update for windkessel model.
    inline void
    windkessel_be_update(double& P_Wk, double& P_l_atrium, const double& Q_l_atrium, const double& t, const double& dt)
    {
        if (t < t_ramp)
            P_Wk       = t * P_ramp / t_ramp; 
        else
            P_Wk       = ((C / dt) * P_Wk + Q_l_atrium) / (C / dt + 1.0 / R_P);
        
        P_l_atrium = P_Wk + R_C * Q_l_atrium;
        return;
    } // windkessel_be_update

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel::CirculationModel(const string& object_name, Pointer<Database> input_db, bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_time(0.0),
      d_nsrc(1),           // number of sets of variables
      d_qsrc(d_nsrc, 0.0), // flux
      d_psrc(d_nsrc, 0.0), // pressure
      d_srcname(d_nsrc),
      d_P_Wk(0.0),
      d_bdry_interface_level_number(numeric_limits<int>::max())
{
#if !defined(NDEBUG)
    assert(!object_name.empty());
#endif
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    if (input_db)
    {
        /*d_bdry_interface_level_number =
            input_db->getIntegerWithDefault("bdry_interface_level_number", d_bdry_interface_level_number);
        string data_time_filename = input_db->getString("data_time_filename");
        string pump_pressure_filename = input_db->getString("pump_pressure_filename");
        string upstream_pressure_filename = input_db->getString("upstream_pressure_filename");
        string downstream_pressure_filename = input_db->getString("downstream_pressure_filename");
        ifstream time_ifs, pump_ifs, inflow_ifs, outflow_ifs;
        time_ifs.open(data_time_filename.c_str(), ios::in);
        while (time_ifs && !time_ifs.eof())
        {
            double t;
            time_ifs >> t;
            data_time.push_back(t);
        }
        pump_ifs.open(upstream_pressure_filename.c_str(), ios::in);
        while (pump_ifs && !pump_ifs.eof())
        {
            double p;
            pump_ifs >> p;
            pump_pressure.push_back(p);
        }
        inflow_ifs.open(upstream_pressure_filename.c_str(), ios::in);
        while (inflow_ifs && !inflow_ifs.eof())
        {
            double p;
            inflow_ifs >> p;
            upstream_pressure.push_back(p);
        }
        outflow_ifs.open(downstream_pressure_filename.c_str(), ios::in);
        while (outflow_ifs && !outflow_ifs.eof())
        {
            double p;
            outflow_ifs >> p;
            downstream_pressure.push_back(p);
        }
        */
        R_P = input_db->getDoubleWithDefault("R_P", R_P); // peripheral resistance (mmHg ml^-1 s)
        R_C = input_db->getDoubleWithDefault("R_C", R_C); // characteristic resistance (mmHg ml^-1 s)
        C   = input_db->getDoubleWithDefault("C",   C);   // total arterial compliance (ml mmHg^-1)
        
        std::cout << "input db got values R_P = " << R_P << "R_C = " << R_C << "C = " << C << "\n";
        
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
        // The downstream (Atrial) pressure is determined by a three-element Windkessel model.
        const double t = d_time;
        double P_l_atrium;
        const double Q_l_atrium = d_qsrc[0];
        double& P_Wk = d_P_Wk;
        windkessel_be_update(P_Wk, P_l_atrium, Q_l_atrium, t, dt);
        d_psrc[0] = P_l_atrium * MMHG_TO_CGS;
    #else 
        // Downstream pressure is fixed to zero 
        d_psrc[0] = 0.0;
    #endif

    // Update the current time.
    d_time += dt;

    // Output the updated values.
    const long precision = plog.precision();
    plog.unsetf(ios_base::showpos);
    plog.unsetf(ios_base::scientific);

    plog << "============================================================================\n"
         << "Circulation model variables at time " << d_time << ":\n";

    plog << "Q_mi\t= ";
    plog.setf(ios_base::showpos);
    plog.setf(ios_base::scientific);
    plog.precision(12);
    plog << -d_qsrc[0] << "ml/s";
    plog << "\n";

    plog << "P_l_atrium\t= ";
    plog.setf(ios_base::showpos);
    plog.setf(ios_base::scientific);
    plog.precision(12);
    plog << d_psrc[0] / MMHG_TO_CGS << " mmHg";
    plog << "\n";

    #ifdef USE_WINDKESSEL
        plog << "P_Wk\t= ";
        plog.setf(ios_base::showpos);
        plog.setf(ios_base::scientific);
        plog.precision(12);
        plog << P_Wk << " mmHg";  // Note that windkessel uses units of mmHg so this needs no conversion
        plog << "\n";
    #endif
    
    plog << "============================================================================\n";

    plog.unsetf(ios_base::showpos);
    plog.unsetf(ios_base::scientific);
    plog.precision(precision);

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
    db->putDouble("d_P_Wk", d_P_Wk);
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
            fout << "% time       "
                 << " P_LA (mmHg) "
                 << " Q_mi (ml/min)"
                 << " P_Wk (mmHg) "
                 << "\n"
                 << "bc_data = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        for (int n = 0; n < d_nsrc; ++n)
        {
            fout << d_time;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(5);
            fout << " " << d_psrc[n] / MMHG_TO_CGS;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(5);
            fout << " " << -d_qsrc[n];
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(5);
            fout << " " << d_P_Wk;
            fout << "; \n";
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
    d_P_Wk = db->getDouble("d_P_Wk");
    d_bdry_interface_level_number = db->getInteger("d_bdry_interface_level_number");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////