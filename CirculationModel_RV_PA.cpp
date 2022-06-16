// Filename: CirculationModel_RV_PA.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_RV_PA.h"
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

CirculationModel_RV_PA::CirculationModel_RV_PA(Pointer<Database> input_db, 
                                               const fourier_series_data *fourier_right_ventricle, 
                                               const fourier_series_data *fourier_right_pa, 
                                               const fourier_series_data *fourier_left_pa, 
                                               string right_ventricle_vertices_file_name,
                                               string right_pa_vertices_file_name,
                                               string left_pa_vertices_file_name,
                                               const double  cycle_duration,
                                               const double  t_offset_bcs_unscaled, 
                                               const double  initial_time, 
                                               double P_initial_pa,
                                               bool rcr_bcs_on,
                                               bool resistance_bcs_on,
                                               bool inductor_bcs_on, 
                                               bool variable_resistance)
    : 
      d_object_name("circ_model_rv_pa"),  // constant name here  
      d_registered_for_restart(true),      // always true
      d_fourier_right_ventricle(fourier_right_ventricle), 
      d_fourier_right_pa(fourier_right_pa),       
      d_fourier_left_pa(fourier_left_pa), 
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_current_idx_series(0),
      d_Q_right_ventricle(0.0), 
      d_Q_right_pa(0.0),
      d_Q_left_pa(0.0),
      d_Q_right_ventricle_previous(0.0),
      d_Q_right_pa_previous(0.0),
      d_Q_left_pa_previous(0.0),
      d_time(initial_time), 
      d_right_pa_P(P_initial_pa), 
      d_right_pa_P_Wk(P_initial_pa),
      d_right_pa_P_distal(P_initial_pa),
      d_right_pa_P_distal_previous(P_initial_pa),
      d_left_pa_P(P_initial_pa),
      d_left_pa_P_Wk(P_initial_pa),
      d_left_pa_P_distal(P_initial_pa),
      d_left_pa_P_distal_previous(P_initial_pa),
      d_area_right_ventricle(0.0),
      d_area_right_pa(0.0),
      d_area_left_pa (0.0),
      d_area_initialized(false), 
      d_rcr_bcs_on(rcr_bcs_on),
      d_resistance_bcs_on(resistance_bcs_on),
      d_inductor_bcs_on(inductor_bcs_on),
      d_variable_resistance(variable_resistance)
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
            d_right_pa_R_proximal = input_db->getDouble("right_pa_R_proximal");
            d_right_pa_R_distal   = input_db->getDouble("right_pa_R_distal");
            d_right_pa_C          = input_db->getDouble("right_pa_C");

            d_left_pa_R_proximal  = input_db->getDouble("left_pa_R_proximal");
            d_left_pa_R_distal    = input_db->getDouble("left_pa_R_distal");
            d_left_pa_C           = input_db->getDouble("left_pa_C");

            std::cout << "input db got values:\n";
            std::cout << "right: R_proximal = " << d_right_pa_R_proximal << "\tR_distal = " << d_right_pa_R_distal << "\tC = " << d_right_pa_C << "\n";
            std::cout << "left : R_proximal = " << d_left_pa_R_proximal << "\tR_distal = " << d_left_pa_R_distal << "\tC = " << d_left_pa_C << "\n";
        }
        else {
            TBOX_ERROR("Must provide valid input_db");
        }
    }

    if (d_resistance_bcs_on){
        d_right_pa_resistance = input_db->getDouble("right_pa_R");
        d_left_pa_resistance  = input_db->getDouble("left_pa_R");
        std::cout << "input db got values:\n";
        std::cout << "right: resistance = " << d_right_pa_resistance << "\n";
        std::cout << "left : R_proximal = " << d_left_pa_resistance << "\n";
    }

    if (d_inductor_bcs_on){
        d_right_pa_inductance = input_db->getDouble("right_pa_L");
        d_left_pa_inductance  = input_db->getDouble("left_pa_L");
        std::cout << "input db got values:\n";
        std::cout << "right: inductance = " << d_right_pa_inductance << "\n";
        std::cout << "left : inductance = " << d_left_pa_inductance << "\n";
    }

    if (d_variable_resistance){
        if (!d_resistance_bcs_on){
            TBOX_ERROR("Must have resistance on to use variable resistance"); 
        }

        // set systole and distole start 
        bool systole_set = false; 
        bool diastole_set = false; 

        double dt = d_fourier_right_ventricle->dt; 

        double t_reduced; 
        double t_scaled; 
        double t_scaled_offset;
        unsigned int k_temp; 
        unsigned int current_idx_series; 
        double pressure_diff_left = 0.0; 
        double pressure_diff_right = 0.0; 
        double pressure_diff_left_prev = 0.0; 
        double pressure_diff_right_prev = 0.0; 
        
        for(double time_temp = 0.0; time_temp < d_cycle_duration; time_temp += dt){

            t_reduced = time_temp - d_cycle_duration * floor(time_temp/d_cycle_duration); 

            // fourier series has its own period, scale to that 
            t_scaled = t_reduced * (d_fourier_right_ventricle->L  / d_cycle_duration); 

            // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
            t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

            // Fourier data here
            // index without periodicity 
            k_temp = (unsigned int) floor(t_scaled_offset / (d_fourier_right_ventricle->dt));
            
            // // take periodic reduction
            current_idx_series = k_temp % (d_fourier_right_ventricle->N_times);

            pressure_diff_left =  d_fourier_right_ventricle->values[current_idx_series] - 
                                          d_fourier_left_pa->values[current_idx_series]; 

            pressure_diff_right = d_fourier_right_ventricle->values[current_idx_series] - 
                                         d_fourier_right_pa->values[current_idx_series]; 


            // negative to positive change 
            if (!systole_set){
                if (((pressure_diff_left_prev < 0.0) && (0.0 <= pressure_diff_left)) || 
                    ((pressure_diff_right_prev < 0.0) && (0.0 <= pressure_diff_right))){
                    d_systole_start = t_reduced; 
                    pout << "found d_systole_start = " << d_systole_start << "\n"; 
                    systole_set = true; 
                }
            }

            if (systole_set && (!diastole_set)){
                if (((pressure_diff_left_prev > 0.0) && (0.0 >= pressure_diff_left)) || 
                    ((pressure_diff_right_prev > 0.0) && (0.0 >= pressure_diff_right))){
                    d_diastole_start = t_reduced; 
                    pout << "found d_diastole_start = " << d_diastole_start << "\n"; 
                    diastole_set = true; 
                }
            }

            pressure_diff_left_prev = pressure_diff_left; 
            pressure_diff_right_prev = pressure_diff_right_prev; 

        }

        if ((!systole_set) || (!diastole_set)){
            TBOX_ERROR("must set systole diastole time"); 
        }

        pout << "Constructor found d_systole_start = " << d_systole_start << ", d_diastole_start = " << d_diastole_start << "\n"; 


    }
    else{
        // not used 
        d_systole_start = 0.0; 
        d_diastole_start = 0.0; 
    }

    if ((d_rcr_bcs_on && d_resistance_bcs_on) || 
        (d_rcr_bcs_on && d_inductor_bcs_on)   ||
        (d_inductor_bcs_on && d_resistance_bcs_on)) {
        TBOX_ERROR("Cannot us two types of bc simulataneously"); 
    }

    double x,x_prev,y,y_prev,z,z_prev; 
    double tol = 1.0e-2; 

    // read vertices from file 
    ifstream right_ventricle_file(right_ventricle_vertices_file_name.c_str(), ios::in);

    if(!right_ventricle_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    right_ventricle_file >> d_n_pts_right_ventricle; 
    
    d_right_ventricle_points_idx1 = new double[d_n_pts_right_ventricle]; 
    d_right_ventricle_points_idx2 = new double[d_n_pts_right_ventricle]; 

    for (int i=0; i<d_n_pts_right_ventricle; i++){
        right_ventricle_file >> x; 
        right_ventricle_file >> d_right_ventricle_points_idx1[i]; 
        right_ventricle_file >> d_right_ventricle_points_idx2[i];
        
        if (i>0){
            if (fabs(x_prev - x) > tol){
                TBOX_ERROR("x coordinates must be consistent\n"); 
            }
        }
        x_prev = x; 

    }
    pout << "to right_ventricle file close\n"; 
    right_ventricle_file.close(); 
    d_right_ventricle_axis = 0; 
    d_right_ventricle_side = 0; 

    // read vertices from file 
    ifstream right_pa_file(right_pa_vertices_file_name.c_str(), ios::in);

    if(!right_pa_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    right_pa_file >> d_n_pts_right_pa; 
    
    d_right_pa_points_idx1 = new double[d_n_pts_right_pa]; 
    d_right_pa_points_idx2 = new double[d_n_pts_right_pa]; 

    for (int i=0; i<d_n_pts_right_pa; i++){
        right_pa_file >> d_right_pa_points_idx1[i]; 
        right_pa_file >> d_right_pa_points_idx2[i]; 
        right_pa_file >> z;
        

        if (i>0){
            if (fabs(z_prev - z) > tol){
                TBOX_ERROR("z coordinates must be consistent\n"); 
            }
        }
        z_prev = z; 

    }
    pout << "to right_pa file close\n"; 
    right_pa_file.close(); 
    d_right_pa_axis = 2; 
    d_right_pa_side = 0; 

    // read vertices from file 
    ifstream left_pa_file(left_pa_vertices_file_name.c_str(), ios::in);

    if(!left_pa_file){
        TBOX_ERROR("Left PA file not found\n"); 
    }

    left_pa_file >> d_n_pts_left_pa; 
    
    d_left_pa_points_idx1 = new double[d_n_pts_left_pa]; 
    d_left_pa_points_idx2 = new double[d_n_pts_left_pa]; 

    for (int i=0; i<d_n_pts_left_pa; i++){
        left_pa_file >> d_left_pa_points_idx1[i]; 
        left_pa_file >> y; 
        left_pa_file >> d_left_pa_points_idx2[i]; 

        if (i>0){
            if (fabs(y_prev - y) > tol){
                TBOX_ERROR("y coordinates must be consistent\n"); 
            }
        }
        y_prev = y; 

    }
    pout << "to left_pa file close\n"; 
    left_pa_file.close();
    d_left_pa_axis = 1; 
    d_left_pa_side = 0; 

    pout << "passed contstructor\n"; 

    return;
} // CirculationModel

CirculationModel_RV_PA::~CirculationModel_RV_PA()
{
    return;
} // ~CirculationModel_RV_PA


void CirculationModel_RV_PA::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    
    double Q_right_ventricle_local = 0.0; 
    double Q_right_pa_local = 0.0; 
    double Q_left_pa_local = 0.0; 

    double area_right_ventricle_local = 0.0; 
    double area_right_pa_local = 0.0; 
    double area_left_pa_local = 0.0; 


    // save old values of Q for taking time derivatives  
    d_Q_right_pa_previous = d_Q_right_pa;
    d_Q_left_pa_previous = d_Q_left_pa;
    d_Q_right_ventricle_previous = d_Q_right_ventricle; 


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

                                const int in_right_ventricle  = this->point_in_right_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_right_pa         = this->point_in_right_pa       (X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_left_pa          = this->point_in_left_pa        (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_right_ventricle && in_right_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }
                                if (in_right_ventricle && in_left_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }
                                if (in_right_pa && in_left_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }

                                if (in_right_ventricle)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_right_ventricle_local += (*U_data)(i_s)* n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_right_ventricle_local += dA;
                                        }

                                    }
                                }

                                if (in_right_pa)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_right_pa_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_right_pa_local += dA;
                                        }

                                    }
                                }

                                if (in_left_pa)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_left_pa_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_left_pa_local += dA;
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

    d_Q_right_ventricle = SAMRAI_MPI::sumReduction(Q_right_ventricle_local);
    d_Q_right_pa        = SAMRAI_MPI::sumReduction(Q_right_pa_local);
    d_Q_left_pa         = SAMRAI_MPI::sumReduction(Q_left_pa_local);

    if (!d_area_initialized){
        d_area_right_ventricle = SAMRAI_MPI::sumReduction(area_right_ventricle_local);
        d_area_right_pa        = SAMRAI_MPI::sumReduction(area_right_pa_local);  
        d_area_left_pa         = SAMRAI_MPI::sumReduction(area_left_pa_local);  
        d_area_initialized = true;       
    }

    // print_summary();

    // bool debug_out_areas = false; 
    // if (debug_out_areas){
    //     pout << "d_area_right_ventricle = " << d_area_right_ventricle << "\n"; 
    //     pout << "d_area_right_pa = " << d_area_right_pa << "\n"; 
    //     pout << "d_area_left_pa = " << d_area_left_pa << "\n"; 
    // }

    d_time += dt; 

    // compute which index in the Fourier series we need here 
    // always use a time in current cycle 
    double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration); 

    // fourier series has its own period, scale to that 
    double t_scaled = t_reduced * (d_fourier_right_ventricle->L  / d_cycle_duration); 

    // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
    double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

    // Fourier data here
    // index without periodicity 
    unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_right_ventricle->dt));
    
    // // take periodic reduction
    d_current_idx_series = k % (d_fourier_right_ventricle->N_times);


    if (d_rcr_bcs_on){
        // The downstream pressure is determined by a three-element Windkessel model.

        double coeff_left = (d_left_pa_C / dt + 1.0 / d_left_pa_R_distal); 
        double coeff_right = (d_right_pa_C / dt + 1.0 / d_right_pa_R_distal);

        // grab the downstream pressures 
        d_right_pa_P_distal_previous = d_right_pa_P_distal; 
        d_right_pa_P_distal = MMHG_TO_CGS * d_fourier_right_pa->values[d_current_idx_series]; 

        d_left_pa_P_distal_previous = d_left_pa_P_distal; 
        d_left_pa_P_distal = MMHG_TO_CGS * d_fourier_left_pa->values[d_current_idx_series]; 

        // hooked to ground version 
        // d_right_pa_P_Wk = ((d_right_pa_C / dt) * d_right_pa_P_Wk + d_Q_right_pa) / (d_right_pa_C / dt + 1.0 / d_right_pa_R_distal);        
        // d_right_pa_P = d_right_pa_P_Wk + d_right_pa_R_proximal * d_Q_right_pa;

        // d_left_pa_P_Wk = ((d_left_pa_C / dt) * d_left_pa_P_Wk + d_Q_left_pa) / (d_left_pa_C / dt + 1.0 / d_left_pa_R_distal);        
        // d_left_pa_P = d_left_pa_P_Wk + d_left_pa_R_proximal * d_Q_left_pa;

        d_right_pa_P_Wk = ((d_right_pa_C / dt) * (d_right_pa_P_Wk - d_right_pa_P_distal_previous) + coeff_right*d_right_pa_P_distal + d_Q_right_pa) / coeff_right;        
        d_right_pa_P = d_right_pa_P_Wk + d_right_pa_R_proximal * d_Q_right_pa;

        d_left_pa_P_Wk = ((d_left_pa_C / dt) * (d_left_pa_P_Wk - d_left_pa_P_distal_previous) + coeff_right*d_left_pa_P_distal + d_Q_left_pa) / coeff_left;        
        d_left_pa_P = d_left_pa_P_Wk + d_left_pa_R_proximal * d_Q_left_pa;

    }
    else if (d_resistance_bcs_on){

        double variable_resistance_coeff = 1.0; 

        if (d_variable_resistance){

            // time to turn resistor on and off 

            // fast off during systole 
            double off_duration_systole  = 0.01;
            // slower on in diastole 
            double on_duration_diastole = 0.1; 

            if (t_reduced < d_systole_start){
                // starts in diastole with closed valve 
                variable_resistance_coeff = 1.0; 
            }
            else if ((d_systole_start <= t_reduced) && (t_reduced < (d_systole_start + off_duration_systole))){
                // ramps up into systole 
                variable_resistance_coeff = 0.5 * (cos( (M_PI/off_duration_systole) * (t_reduced - d_systole_start)  ) + 1); 
            }
            else if (((d_systole_start + off_duration_systole) <= t_reduced) && (t_reduced < d_diastole_start)){
                // off during systole 
                variable_resistance_coeff = 0.0; 
            }
            else if ((d_diastole_start <= t_reduced) && (t_reduced < (d_diastole_start + on_duration_diastole))){
                // ramps on in diastole 
                variable_resistance_coeff = 0.5 * (-cos( (M_PI/on_duration_diastole) * (t_reduced - d_diastole_start)  ) + 1); 
            }
            else {
                // then on for remainder of cycle 
                variable_resistance_coeff = 1.0; 
            }
        
        }


        // pressure upstream of resistance determined by series 
        d_right_pa_P_Wk = MMHG_TO_CGS * d_fourier_right_pa->values[d_current_idx_series]; 
        d_left_pa_P_Wk  = MMHG_TO_CGS * d_fourier_left_pa->values[d_current_idx_series]; 

        // resistance bcs determine outlet pressure 
        d_right_pa_P = d_right_pa_P_Wk + d_right_pa_resistance * variable_resistance_coeff * d_Q_right_pa;
        d_left_pa_P  = d_left_pa_P_Wk  + d_left_pa_resistance  * variable_resistance_coeff * d_Q_left_pa;

    }

    else if (d_inductor_bcs_on){
        // pressure upstream of resistance determined by series 
        d_right_pa_P_Wk = MMHG_TO_CGS * d_fourier_right_pa->values[d_current_idx_series]; 
        d_left_pa_P_Wk  = MMHG_TO_CGS * d_fourier_left_pa->values[d_current_idx_series]; 

        // inductance bcs determine update on pressure 
        d_right_pa_P = d_right_pa_P_Wk + d_right_pa_inductance * (d_Q_right_pa - d_Q_right_pa_previous)/dt;
        d_left_pa_P  = d_left_pa_P_Wk  + d_left_pa_inductance  * (d_Q_left_pa  - d_Q_left_pa_previous )/dt ;
    }

    else {
        d_right_pa_P = MMHG_TO_CGS * d_fourier_right_pa->values[d_current_idx_series]; 
        d_left_pa_P  = MMHG_TO_CGS * d_fourier_left_pa->values[d_current_idx_series]; 
    }

    // bool debug_out = false; 
    // if (debug_out){
    //     pout << "circ mode: d_time = " << d_time << ", d_current_idx_series = " << d_current_idx_series << "\n"; 
    //     pout << "t_reduced = " << t_reduced << " t_scaled = " << t_scaled << " t_scaled_offset = " << t_scaled_offset << "\n"; 
    //     pout << "k (unreduced idx) = " << k << " d_current_idx_series = " << d_current_idx_series << "\n\n"; 
    // }


    writeDataFile(); 

} // advanceTimeDependentData

void CirculationModel_RV_PA::set_Q_valve(double Q_valve){
    d_Q_valve = Q_valve; 
}



void
CirculationModel_RV_PA::putToDatabase(Pointer<Database> db)
{

    db->putInteger("d_current_idx_series", d_current_idx_series); 
    db->putDouble("d_Q_right_ventricle", d_Q_right_ventricle); 
    db->putDouble("d_Q_right_pa", d_Q_right_pa);
    db->putDouble("d_Q_left_pa", d_Q_left_pa);
    db->putDouble("d_Q_right_ventricle_previous", d_Q_right_ventricle_previous);
    db->putDouble("d_Q_right_pa_previous", d_Q_right_pa_previous);
    db->putDouble("d_Q_left_pa_previous", d_Q_left_pa_previous);
    db->putDouble("d_Q_valve", d_Q_valve);
    db->putDouble("d_right_pa_P", d_right_pa_P);
    db->putDouble("d_right_pa_P_Wk", d_right_pa_P_Wk);
    db->putDouble("d_right_pa_P_distal", d_right_pa_P_distal);
    db->putDouble("d_right_pa_P_distal_previous", d_right_pa_P_distal_previous);
    db->putDouble("d_left_pa_P",d_left_pa_P);
    db->putDouble("d_left_pa_P_Wk",d_left_pa_P_Wk);
    db->putDouble("d_left_pa_P_distal", d_left_pa_P_distal);
    db->putDouble("d_left_pa_P_distal_previous", d_left_pa_P_distal_previous);
    db->putDouble("d_time", d_time); 
    db->putBool("d_rcr_bcs_on", d_rcr_bcs_on); 
    db->putBool("d_resistance_bcs_on", d_resistance_bcs_on); 
    db->putBool("d_inductor_bcs_on", d_inductor_bcs_on); 
    db->putBool("d_variable_resistance", d_variable_resistance); 
    return; 
} // putToDatabase

void CirculationModel_RV_PA::print_summary(){

    double P_right_ventricle = d_fourier_right_ventricle->values[d_current_idx_series]; 
    double P_right_pa; 
    double P_left_pa; 

    if (d_rcr_bcs_on){
        P_right_pa        = d_right_pa_P / MMHG_TO_CGS;
        P_left_pa         = d_left_pa_P / MMHG_TO_CGS;
    }
    else{
        P_right_pa        = d_fourier_right_pa->values[d_current_idx_series];
        P_left_pa         = d_fourier_left_pa->values[d_current_idx_series];
    }

    pout << "rcr_bcs_on = " << d_rcr_bcs_on << "\n"; 
    pout << "% time \t P_right_ventricle (mmHg)\t P_right_pa (mmHg)\t P_left_pa (mmHg)\t Q_right_ventricle (ml/s)\t d_Q_right_pa (ml/s)\t d_Q_right_pa (ml/s)\tQ_valve (ml/s) \t idx" ;
    if (d_rcr_bcs_on || d_resistance_bcs_on){
        pout << "\t right_pa_P_Wk\t left_pa_P_Wk\t "; 
    }
    pout << "\n";
    pout << d_time << " " << P_right_ventricle <<  " " << P_right_pa << " " << P_left_pa << " " << d_Q_right_ventricle << " " << d_Q_right_pa << " " << d_Q_left_pa << " " << d_Q_valve << " " << d_current_idx_series; 
    if (d_rcr_bcs_on || d_resistance_bcs_on){
        pout  << " " << d_right_pa_P_Wk << " " << d_left_pa_P_Wk; 
    }
    pout << "\n";

}

int CirculationModel_RV_PA::point_in_right_ventricle(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_right_ventricle_axis) || (side != d_right_ventricle_side))
        return 0; 

    return pnpoly(d_n_pts_right_ventricle, d_right_ventricle_points_idx1, d_right_ventricle_points_idx2, testx, testy); 
}

int CirculationModel_RV_PA::point_in_right_pa(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_right_pa_axis) || (side != d_right_pa_side))
        return 0; 

    return pnpoly(d_n_pts_right_pa, d_right_pa_points_idx1, d_right_pa_points_idx2, testx, testy); 
}

int CirculationModel_RV_PA::point_in_left_pa(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_left_pa_axis) || (side != d_left_pa_side))
        return 0; 

    return pnpoly(d_n_pts_left_pa, d_left_pa_points_idx1, d_left_pa_points_idx2, testx, testy); 
}

void CirculationModel_RV_PA::write_plot_code()
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
        fout << "times   =  bc_vals(:,1);\n";
        fout << "p_rv    =  bc_vals(:,2);\n";
        fout << "p_rpa   =  bc_vals(:,3);\n";
        fout << "p_lpa   =  bc_vals(:,4);\n";
        fout << "q_rv    = -bc_vals(:,5); \n";
        fout << "q_rpa   =  bc_vals(:,6);\n";
        fout << "q_lpa   =  bc_vals(:,7);\n";
        fout << "q_valve =  bc_vals(:,8);\n";
        fout << "p_wk_rpa =  bc_vals(:,9);\n";
        fout << "p_wk_lpa =  bc_vals(:,10);\n";
        fout << "load '../bc_variables_experimental.mat'\n";
        fout << "subplot(2,1,1)\n";
        fout << "plot(times, p_rv, 'k')\n";
        fout << "hold on\n";
        fout << "plot(times, p_rpa, ':k')\n";
        fout << "plot(times, p_lpa, '-.k')\n";
        fout << "plot(times_exp, p_rv_exp)\n";
        fout << "plot(times_exp, p_pa_exp)\n";
        fout << "plot(times, p_wk_rpa)\n";
        fout << "plot(times, p_wk_lpa)\n";
        fout << "% legend('P RV', 'P RPA', 'PWK RPA', 'P LPA', 'PWK LPA', 'P EXP RV', 'P EXP PA', Location','NorthEastOutside');\n";
        fout << "legend('P RV', 'P RPA', 'P LPA', 'P EXP RV', 'P EXP PA', 'WK RPA', 'WK LPA', 'Location','NorthEastOutside');\n";
        fout << "xlabel('t (s)')\n";
        fout << "ylabel('P (mmHg)')\n";
        fout << "subplot(2,1,2)\n";
        fout << "plot(times, q_rv, 'k')\n";
        fout << "hold on\n";
        fout << "plot(times, q_rpa, '--k')\n";
        fout << "plot(times, q_lpa, '-.k')\n";
        fout << "plot(times_two_cycles, q_rpa_exp)\n";
        fout << "plot(times_two_cycles, q_lpa_exp)\n";
        fout << "plot(times_two_cycles, q_rv_exp)\n";
        fout << "plot(bc_vals(:,1), zeros(size(q_rv)), ':k')\n";
        fout << "legend('Q RV', 'Q RPA', 'Q LPA', 'Q EXP RV', 'Q EXP RPA', 'Q EXP LPA', 'Location', 'NorthEastOutside')\n";
        fout << "xlabel('t (s)')\n";
        fout << "ylabel('Flow (ml/s)')\n";
        fout << "set(fig, 'Position', [100, 100, 1000, 750])\n";
        fout << "set(fig,'PaperPositionMode','auto')\n";
        fout << "printfig(fig, 'bc_model_variables_with_experimental')\n";
        fout << "q_mean = mean(q_rv)\n";
    }
    return;
}



/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
    CirculationModel_RV_PA::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time \t P_right_ventricle (mmHg)\t P_right_pa (mmHg)\t P_left_pa (mmHg)\t d_Q_right_ventricle (ml/s)\t d_Q_right_pa (ml/s)\td_Q_left_pa (ml/s) \td_Q_valve (ml/s) \t d_right_pa_P_Wk \t d_left_pa_P_Wk"; 
            if (d_rcr_bcs_on){
                fout << "d_right_pa_P_distal \t d_left_pa_P_distal"; 
            }
            fout << "\n"
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);

        fout << d_time;
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        double P_right_ventricle = d_fourier_right_ventricle->values[d_current_idx_series]; 

        fout << " " << P_right_ventricle <<  " " << d_right_pa_P/MMHG_TO_CGS << " " << d_left_pa_P/MMHG_TO_CGS;
        fout << " " << d_Q_right_ventricle << " " << d_Q_right_pa << " " << d_Q_left_pa << " " << d_Q_valve;         
        fout << " " << d_right_pa_P_Wk/MMHG_TO_CGS << " " << d_left_pa_P_Wk/MMHG_TO_CGS;
        if(d_rcr_bcs_on){
            fout << " " << d_right_pa_P_distal/MMHG_TO_CGS << " " << d_left_pa_P_distal/MMHG_TO_CGS;
        }
        fout << "; \n";

    }

    return;
} // writeDataFile

void
CirculationModel_RV_PA::getFromRestart()
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

    d_current_idx_series         = db->getInteger("d_current_idx_series"); 
    d_Q_right_ventricle          = db->getDouble("d_Q_right_ventricle"); 
    d_Q_right_pa                 = db->getDouble("d_Q_right_pa");
    d_Q_left_pa                  = db->getDouble("d_Q_left_pa");
    d_Q_right_ventricle          = db->getDouble("d_Q_right_ventricle_previous"); 
    d_Q_right_pa_previous        = db->getDouble("d_Q_right_pa_previous");
    d_Q_left_pa_previous         = db->getDouble("d_Q_left_pa_previous");
    d_Q_valve                    = db->getDouble("d_Q_valve");
    d_right_pa_P                 = db->getDouble("d_right_pa_P");
    d_right_pa_P_Wk              = db->getDouble("d_right_pa_P_Wk");
    d_right_pa_P_distal          = db->getDouble("d_right_pa_P_distal");
    d_right_pa_P_distal_previous = db->getDouble("d_right_pa_P_distal_previous");
    d_left_pa_P                  = db->getDouble("d_left_pa_P");
    d_left_pa_P_Wk               = db->getDouble("d_left_pa_P_Wk");
    d_left_pa_P_distal           = db->getDouble("d_left_pa_P_distal");
    d_left_pa_P_distal_previous  = db->getDouble("d_left_pa_P_distal_previous");
    d_time                       = db->getDouble("d_time");
    d_rcr_bcs_on                 = db->getBool("d_rcr_bcs_on"); 
    d_resistance_bcs_on          = db->getBool("d_resistance_bcs_on"); 
    d_inductor_bcs_on            = db->getBool("d_inductor_bcs_on"); 
    d_variable_resistance        = db->getBool("d_variable_resistance"); 
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////