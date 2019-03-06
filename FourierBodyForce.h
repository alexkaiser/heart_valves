
// FourierBodyForce.cpp
// Created Alexander D. Kaiser, 6/2016

// Modified from:
// Filename: FeedbackForcer.h
// Created on 08 Sep 2007 by Boyce Griffith

#ifndef included_FourierBodyForce
#define included_FourierBodyForce

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES/
#include <petscsys.h>

#include <boundary_condition_util.h>

// IBAMR INCLUDES
#include <ibamr/INSHierarchyIntegrator.h>

// NAMESPACE
#include <ibamr/app_namespaces.h>


/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FourierBodyForce.cpp is an implementation of the strategy class
 * CartGridFunction that is used to specify a time dependent body force in the z direction
 *
 * The forces is given -e_z * (delta P) / Length 
 * where delta P is the pressure drop, Length is the domain size
 */
class FourierBodyForce : public CartGridFunction
{
 public:
  /*!
   * \brief Constructor
   */
  FourierBodyForce(const fourier_series_data* fourier,
                   const bool use_circ_model,
                   CirculationModel *circ_model,
		           const INSHierarchyIntegrator* fluid_solver,
		           Pointer<PatchHierarchy<NDIM> > patch_hierarchy);

  /*!
   * \brief Destructor.
   */
  virtual ~FourierBodyForce();

  /*!
   * \name Implementation of CartGridFunction interface.
   */
  //\{

  /*!
   * \brief Indicates whether the concrete CartGridFunction object is
   * time-dependent.
   */
  bool isTimeDependent() const;

  /*!
   * \brief Set data on the specified patch interior.
   */
  void setDataOnPatch(int data_idx,
		      Pointer<hier::Variable<NDIM> > var,
		      Pointer<Patch<NDIM> > patch,
		      double data_time,
		      bool initial_time = false,
		      Pointer<PatchLevel<NDIM> > patch_level = Pointer<PatchLevel<NDIM> >(NULL));

  //\}

  const fourier_series_data *d_fourier;
  const bool d_use_circ_model; 
  CirculationModel *d_circ_model;
  double d_flux_z; 

 private:
  /*!
   * \brief Default constructor.
   *
   * \note This constructor is not implemented and should not be used.
   */
  FourierBodyForce();

  /*!
   * \brief Copy constructor.
   *
   * \note This constructor is not implemented and should not be used.
   *
   * \param from The value to copy to this object.
   */
  FourierBodyForce(const FourierBodyForce& from);

  /*!
   * \brief Assignment operator.
   *
   * \note This operator is not implemented and should not be used.
   *
   * \param that The value to assign to this object.
   *
   * \return A reference to this object.
   */
  FourierBodyForce& operator=(const FourierBodyForce& that);

  const INSHierarchyIntegrator* const d_fluid_solver;
  Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
  
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <FourierBodyForce.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FourierBodyForce
