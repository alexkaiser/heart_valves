// Filename: VelocityBcCoefs.h
// Created on 04 May 2007 by Boyce Griffith

// Modified for Fourier series input data 
// 5/2016, Alex Kaiser 

#ifndef included_VelocityBcCoefs
#define included_VelocityBcCoefs

// PETSC INCLUDES
#include <petscsys.h>

// SAMRAI INCLUDES
#include <RobinBcCoefStrategy.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>

#define MMHG_TO_CGS 1333.22368

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
    VelocityBcCoefs(const fourier_series_data *fourier);

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


















