// Filename: VelocityBcCoefs.h
// Created on 04 May 2007 by Boyce Griffith

// Modified for Fourier series input data 
// 5/2016, Alexander D. Kaiser 



// PETSC INCLUDES
#include <petscsys.h>

// SAMRAI INCLUDES
#include <RobinBcCoefStrategy.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>

// #include "CirculationModel_with_lv.h"
// forward declare this class rather than include to remove circular references 
class CirculationModel_with_lv; 
class CirculationModel_RV_PA; 

#include "CirculationModel.h"



#define MMHG_TO_CGS 1333.22368


#ifndef included_fourier_series_data
#define included_fourier_series_data

class fourier_series_data{
    
    public: 
    
        unsigned int N_times; 
        double L; 
        double *values; 
        double dt; 

        // Only implemented constructor 
        fourier_series_data(string file_name, double dt_input); 
        
        // Destructor 
        ~fourier_series_data(); 
        
        void print_values(); 

    private: 

}; 

#endif 


#ifndef included_VelocityBcCoefs
#define included_VelocityBcCoefs

/*!
 * \brief Class VelocityBcCoefs is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs : public RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs(const fourier_series_data *fourier, CirculationModel *circ_model);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                    Pointer<ArrayData<NDIM, double> >& bcoef_data,
                    Pointer<ArrayData<NDIM, double> >& gcoef_data,
                    const Pointer<hier::Variable<NDIM> >& variable,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
     IntVector<NDIM> numberOfExtensionsFillable() const;

    const fourier_series_data *d_fourier;
    CirculationModel *d_circ_model;

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    // VelocityBcCoefs(const VelocityBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs& operator=(const VelocityBcCoefs& that);

};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <VelocityBcCoefs.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityBcCoefs





#ifndef included_VelocityBcCoefs_lv_aorta
#define included_VelocityBcCoefs_lv_aorta

/*!
 * \brief Class VelocityBcCoefs_lv_aorta is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs_lv_aorta : public RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs_lv_aorta(const int comp_idx,
                             CirculationModel_with_lv* circ_model_with_lv);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs_lv_aorta();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                    Pointer<ArrayData<NDIM, double> >& bcoef_data,
                    Pointer<ArrayData<NDIM, double> >& gcoef_data,
                    const Pointer<hier::Variable<NDIM> >& variable,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
     IntVector<NDIM> numberOfExtensionsFillable() const;

    const int d_comp_idx; 
    CirculationModel_with_lv* d_circ_model_with_lv; 

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    // VelocityBcCoefs_lv_aorta(const VelocityBcCoefs_lv_aorta& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs_lv_aorta& operator=(const VelocityBcCoefs_lv_aorta& that);

};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <VelocityBcCoefs_lv_aorta.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityBcCoefs_lv_aorta





#ifndef included_VelocityBcCoefs_RV_PA
#define included_VelocityBcCoefs_RV_PA

/*!
 * \brief Class VelocityBcCoefs_RV_PA is an implementation of the strategy class
 * RobinBcCoefStrategy that is used to specify velocity boundary conditions that
 * are determined by a circulation model.
 */
class VelocityBcCoefs_RV_PA : public RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    VelocityBcCoefs_RV_PA(const int comp_idx,
                          CirculationModel_RV_PA* circ_model_rv_pa);

    /*!
     * \brief Destructor.
     */
    virtual ~VelocityBcCoefs_RV_PA();

    /*!
     * \name Implementation of RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     */
    void setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                    Pointer<ArrayData<NDIM, double> >& bcoef_data,
                    Pointer<ArrayData<NDIM, double> >& gcoef_data,
                    const Pointer<hier::Variable<NDIM> >& variable,
                    const Patch<NDIM>& patch,
                    const BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     */
     IntVector<NDIM> numberOfExtensionsFillable() const;

    const int d_comp_idx; 
    CirculationModel_RV_PA* d_circ_model_rv_pa; 

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    // VelocityBcCoefs_RV_PA(const VelocityBcCoefs_RV_PA& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VelocityBcCoefs_RV_PA& operator=(const VelocityBcCoefs_RV_PA& that);

};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <VelocityBcCoefs_RV_PA.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VelocityBcCoefs_RV_PA
















