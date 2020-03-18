// Filename: CirculationModel.h
// Created on 04 May 2007 by Boyce Griffith

#ifndef included_CirculationModel
#define included_CirculationModel

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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class CirculationModel
 */
class CirculationModel : public Serializable
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
     * \brief Windkessel model data.
     */
    double d_time;
    int d_nsrc;
    vector<double> d_qsrc, d_psrc;
    vector<string> d_srcname;
    double d_P_Wk;
    double d_p_opposite;

    double d_R_proximal; 
    double d_R_distal; 
    double d_C; 


    /*!
     * \brief The level of the patch hierarchy on which the Lagrangian
     * structures that interface the boundary are located.
     */
    int d_bdry_interface_level_number;

    /*!
     * \brief Constructor
     */
    CirculationModel(const string& object_name, Pointer<Database> input_db, bool register_for_restart, double P_initial);

    /*!
     * \brief Destructor.
     */
    virtual ~CirculationModel();

    void windkessel_be_update(double& P_Wk, double& P_boundary, const double& Q_l_atrium, const double& dt); 

    /*!
     * \brief Advance time-dependent data.
     */
    void advanceTimeDependentData(const double dt, const double Q_input);

    /*!
     * \name Implementation of Serializable interface.
     */

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database point must be non-null.
     */
    void putToDatabase(Pointer<Database> db);

    void write_plot_code(); 


private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CirculationModel(const CirculationModel& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CirculationModel& operator=(const CirculationModel& that);

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

#endif //#ifndef included_CirculationModel
