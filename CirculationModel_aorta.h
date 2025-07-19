// Filename: CirculationModel_aorta.h
// Created on 04 May 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser


#ifndef included_CirculationModel_aorta
#define included_CirculationModel_aorta

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <tbox/Database.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <vector>

// NAMESPACE
#include <ibamr/app_namespaces.h>

#include <boundary_condition_util.h>

#include "pnpoly.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel_aorta
 */
class CirculationModel_aorta : public Serializable
{
public:
    /*!
     * \brief The object name.
     */
    string d_object_name;

    /*!
     * \brief Whether the object is registered with the restart manager.
     */
    bool d_registered_for_restart;

    /*!
     * \brief model data.
     */
    const fourier_series_data *d_fourier_ventricle;
    int     d_n_pts_ventricle;
    double* d_ventricle_points_idx1;
    double* d_ventricle_points_idx2;
    int     d_ventricle_axis; 
    int     d_ventricle_side; 
    double  d_ventricle_P;

    ventricle_0D_model *d_ventricle_0D;     

    bool d_rcr_bcs_on; 
    bool d_ventricle_0D_on; 
    bool d_lv_systolic_on;

    // const fourier_series_data *d_fourier_aorta;
    int     d_n_pts_aorta;
    double* d_aorta_points_idx1;
    double* d_aorta_points_idx2;
    int     d_aorta_axis; 
    int     d_aorta_side;
    double  d_aorta_P;
    double  d_aorta_P_Wk;
    double  d_aorta_R_proximal; 
    double  d_aorta_R_distal; 
    double  d_aorta_C;

    double  d_cycle_duration;
    double  d_t_offset_bcs_unscaled;
    unsigned int  d_current_idx_series; 
    double        d_Q_ventricle; 
    double        d_Q_aorta;
    double        d_Q_valve;
    double        d_time; 
    double        d_area_ventricle; 
    double        d_area_aorta; 
    double        d_area_left_pa; 
    bool          d_area_initialized;

    double d_p_extender_mean; 
    double d_p_extender_point; 

    bool d_P_initial_aorta_equal_to_ventricle; 
    double d_p_equal_fraction; 
    double d_P_min_linear_interp; 
    double d_rcr_on_time; 

    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief Constructor
     */
    CirculationModel_aorta(Pointer<Database> input_db, 
                                               const fourier_series_data *fourier_ventricle, 
                                               string ventricle_vertices_file_name,
                                               string aorta_vertices_file_name,
                                               const double  cycle_duration,
                                               const double  t_offset_bcs_unscaled, 
                                               const double  initial_time, 
                                               double P_initial_aorta,
                                               bool rcr_bcs_on,
                                               bool ventricle_0D_on, 
                                               bool P_initial_aorta_equal_to_ventricle,
                                               double rcr_on_time);  


    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel_aorta();

    /*!
     * \brief Advance time-dependent data.
     */
    void advanceTimeDependentData(const double dt,
                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int U_idx,
                                  const int P_idx,
                                  const int wgt_cc_idx,
                                  const int wgt_sc_idx);

    void set_Q_valve(double Q_valve); 

    void set_extender_pressures(double p_extender_mean, double p_extender_point); 

    /*!
     * \name Implementation of Serializable interface.
     */

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database point must be non-null.
     */
    void putToDatabase(Pointer<Database> db);

    // basic data summary to stdout 
    void print_summary(); 
    void print_bc_debug();

    int point_in_ventricle(double testx, double testy, int axis, int side);
    int point_in_aorta(double testx, double testy, int axis, int side);

    void write_plot_code(); 

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel_aorta(const CirculationModel_aorta& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel_aorta& operator=(const CirculationModel_aorta& that);

    /*!
     * Write out source/sink state data to disk.
     */
    void writeDataFile() const;

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    void getFromRestart();
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <CirculationModel.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CirculationModel_aorta
