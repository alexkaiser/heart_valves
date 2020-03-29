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

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

#define MMHG_TO_CGS 1333.22368

namespace
{
    // Name of output file.
    static const string DATA_FILE_NAME = "bc_data.m";

}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel::CirculationModel(const string& object_name, Pointer<Database> input_db, bool register_for_restart, double P_initial)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_time(0.0),
      d_nsrc(1),           // number of sets of variables
      d_qsrc(d_nsrc, 0.0), // flux
      d_psrc(d_nsrc, P_initial), // pressure
      d_p_opposite(0.0), // pressure opposite face 
      d_srcname(d_nsrc),
      d_P_Wk(P_initial),
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
        d_R_proximal = input_db->getDouble("R_proximal"); 
        d_R_distal   = input_db->getDouble("R_distal"); 
        d_C          = input_db->getDouble("C");
        
        std::cout << "input db got values R_proximal = " << d_R_proximal << "\tR_distal = " << d_R_distal << "\tC = " << d_C << "\n";   
    }
        else
    {
        TBOX_ERROR("Must provide valid input_db"); 
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
        //           (1) windkessel
        d_srcname[0] = "windkessel       ";        
    }
    return;
} // CirculationModel

CirculationModel::~CirculationModel()
{
    return;
} // ~CirculationModel


void CirculationModel::windkessel_be_update(double& P_Wk, double& P_boundary, const double& Q_l_atrium, const double& dt)
{
    // Backward Euler update for windkessel model.
    P_Wk       = ((d_C / dt) * P_Wk + Q_l_atrium) / (d_C / dt + 1.0 / d_R_distal);        
    P_boundary = P_Wk + d_R_proximal * Q_l_atrium;
    return;
} // windkessel_be_update


void
CirculationModel::advanceTimeDependentData(const double dt, const double Q_input)
{
    d_qsrc[0] = Q_input; 

    // The downstream pressure is determined by a three-element Windkessel model.
    double P_boundary;
    const double Q_source = d_qsrc[0];
    double& P_Wk = d_P_Wk;
    windkessel_be_update(P_Wk, P_boundary, Q_source, dt);
    d_psrc[0] = P_boundary;


    // Update the current time.
    d_time += dt;

    // Output the updated values.
    const long precision = plog.precision();
    plog.unsetf(ios_base::showpos);
    plog.unsetf(ios_base::scientific);

    plog << "============================================================================\n"
         << "Circulation model variables at time " << d_time << ":\n";

    plog.setf(ios_base::showpos);
    plog.setf(ios_base::scientific);
    plog.precision(5);

    plog << "Q           = " << d_qsrc[0] << " ml/s\n";
    plog << "P_boundary  = " << d_psrc[0]/MMHG_TO_CGS << " mmHg\t" << d_psrc[0] << " dynes/cm^2\n";
    plog << "P_Wk        = " << d_psrc[0]/MMHG_TO_CGS << " mmHg\t" << P_Wk      << " dynes/cm^2\n" ;  
    
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
    db->putDouble("d_p_opposite", d_p_opposite);
    db->putInteger("d_bdry_interface_level_number", d_bdry_interface_level_number);
    return;
} // putToDatabase


void CirculationModel::write_plot_code()
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);
        fout << "];\n";  
        fout << "MMHG_TO_CGS = 1333.22368;\n"; 
        fout << "fig = figure;\n";
        fout << "times = bc_vals(:,1);\n"; 
        fout << "p_aorta = bc_vals(:,2)/MMHG_TO_CGS;\n"; 
        fout << "q_aorta = bc_vals(:,3);\n"; 
        fout << "p_wk = bc_vals(:,4)/MMHG_TO_CGS;\n"; 
        fout << "p_lv = bc_vals(:,5)/MMHG_TO_CGS;\n"; 
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
        fout << "plot(bc_vals(:,1), net_flux, '--k')\n";
        fout << "plot(bc_vals(:,1), 0*net_flux, ':k')\n";
        fout << "legend('Q', 'net Q', 'Location','NorthEastOutside')\n";
        fout << "xlabel('t (s)')\n";
        fout << "ylabel('Flow (ml/s), Net Flow (ml)')\n";
        fout << "set(fig, 'Position', [100, 100, 1000, 750])\n";
        fout << "set(fig,'PaperPositionMode','auto')\n";
        fout << "printfig(fig, 'bc_model_variables')\n";
        fout << "mean_p_aorta = mean(p_aorta)\n";
        fout << "mean_p_wk    = mean(p_wk)\n";
        fout << "mean_p_lv    = mean(p_lv)\n";
    }
    return;
}


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
                 << " P_outlet (dynes/cm^2)"
                 << " Q (ml/s)"
                 << " P_Wk (dynes/cm^2)"
                 << " P_opposite (dynes/cm^2)" 
                 << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        static ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        for (int n = 0; n < d_nsrc; ++n)
        {
            fout << d_time;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(10);
            fout << " " << d_psrc[n];
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(10);
            fout << " " << d_qsrc[n];
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(10);
            fout << " " << d_P_Wk;
            fout.setf(ios_base::scientific);
            fout.setf(ios_base::showpos);
            fout.precision(10);
            fout << " " << d_p_opposite;
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
    d_p_opposite = db->getDouble("d_p_opposite");
    d_bdry_interface_level_number = db->getInteger("d_bdry_interface_level_number");
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////