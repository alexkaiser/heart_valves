// Filename: CirculationModel_with_lv.h
// Created on 04 May 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser


#ifndef included_CirculationModel_with_lv
#define included_CirculationModel_with_lv

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


int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy); 



/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel_with_lv
 */
class CirculationModel_with_lv : public Serializable
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
    const fourier_series_data *d_fourier_aorta;
    const fourier_series_data *d_fourier_atrium;
    const fourier_series_data *d_fourier_ventricle; 
    int     d_n_pts_aorta;
    double* d_aorta_points_x;
    double* d_aorta_points_y;
    int     d_n_pts_atrium;
    double* d_atrium_points_x;
    double* d_atrium_points_y;
    double  d_cycle_duration;
    double  d_t_offset_bcs_unscaled;
    unsigned int  d_current_idx_series; 
    double        d_Q_aorta; 
    double        d_Q_left_atrium;
    double        d_Q_mitral;
    double        d_time; 

    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief Constructor
     */
    CirculationModel_with_lv(const fourier_series_data *fourier_aorta, 
                             const fourier_series_data *fourier_atrium,
                             const fourier_series_data *fourier_ventricle,
                             string aorta_vertices_file_name,
                             string atrium_vertices_file_name,
                             const double  cycle_duration,
                             const double  t_offset_bcs_unscaled,
                             const double  initial_time); 

    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel_with_lv();

    /*!
     * \brief Advance time-dependent data.
     */
    void advanceTimeDependentData(const double dt,
                                  const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                  const int U_idx,
                                  const int P_idx,
                                  const int wgt_cc_idx,
                                  const int wgt_sc_idx);

    void set_Q_mitral(double Q_mitral); 


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

    int point_in_aorta(double testx, double testy); 

    int point_in_atrium(double testx, double testy); 


private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel_with_lv(const CirculationModel_with_lv& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel_with_lv& operator=(const CirculationModel_with_lv& that);

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

#endif //#ifndef included_CirculationModel_with_lv
