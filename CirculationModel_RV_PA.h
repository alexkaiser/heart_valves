// Filename: CirculationModel_RV_PA.h
// Created on 04 May 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser


#ifndef included_CirculationModel_RV_PA
#define included_CirculationModel_RV_PA

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
 * \brief Class CirculationModel_RV_PA
 */
class CirculationModel_RV_PA : public Serializable
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
    const fourier_series_data *d_fourier_right_ventricle;
    int     d_n_pts_right_ventricle;
    double* d_right_ventricle_points_idx1;
    double* d_right_ventricle_points_idx2;
    int     d_right_ventricle_axis; 
    int     d_right_ventricle_side; 

    const fourier_series_data *d_fourier_right_pa;
    int     d_n_pts_right_pa;
    double* d_right_pa_points_idx1;
    double* d_right_pa_points_idx2;
    int     d_right_pa_axis; 
    int     d_right_pa_side;

    const fourier_series_data *d_fourier_left_pa; 
    int     d_n_pts_left_pa;
    double* d_left_pa_points_idx1;
    double* d_left_pa_points_idx2;
    int     d_left_pa_axis; 
    int     d_left_pa_side;

    double  d_cycle_duration;
    double  d_t_offset_bcs_unscaled;
    unsigned int  d_current_idx_series; 
    double        d_Q_right_ventricle; 
    double        d_Q_right_pa;
    double        d_Q_left_pa;
    double        d_Q_valve;
    double        d_time; 
    double        d_area_right_ventricle; 
    double        d_area_right_pa; 
    double        d_area_left_pa; 
    bool          d_area_initialized;

    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief Constructor
     */
    CirculationModel_RV_PA(const fourier_series_data *fourier_right_ventricle, 
                                                       const fourier_series_data *fourier_right_pa, 
                                                       const fourier_series_data *fourier_left_pa, 
                                                       string right_ventricle_vertices_file_name,
                                                       string right_pa_vertices_file_name,
                                                       string left_pa_vertices_file_name,
                                                       const double  cycle_duration,
                                                       const double  t_offset_bcs_unscaled, 
                                                       const double  initial_time); 

    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel_RV_PA();

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

    int point_in_right_ventricle(double testx, double testy, int axis, int side);
    int point_in_right_pa(double testx, double testy, int axis, int side);
    int point_in_left_pa(double testx, double testy, int axis, int side);

    void write_plot_code(); 

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel_RV_PA(const CirculationModel_RV_PA& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel_RV_PA& operator=(const CirculationModel_RV_PA& that);

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

#endif //#ifndef included_CirculationModel_RV_PA
