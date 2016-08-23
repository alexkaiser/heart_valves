// Filename: FeedbackForcer.h
// Created on 08 Sep 2007 by Boyce Griffith

#ifndef included_FeedbackForcer
#define included_FeedbackForcer

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

#include <boundary_condition_util.h>

// IBAMR INCLUDES
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>


// NAMESPACE
#include <ibamr/app_namespaces.h>



/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class FeedbackForcer is an implementation of the strategy class
 * CartGridFunction that is used to specify velocity boundary conditions via a
 * feedback forcing (i.e., penalty) method.
 */
class FeedbackForcer : public CartGridFunction
{
 public:
  /*!
   * \brief Constructor
   */
  FeedbackForcer(const VelocityBcCoefs* velocity_bc,
		 const INSHierarchyIntegrator* fluid_solver,
		 Pointer<PatchHierarchy<NDIM> > patch_hierarchy);

  /*!
   * \brief Destructor.
   */
  virtual ~FeedbackForcer();

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



 private:
  /*!
   * \brief Default constructor.
   *
   * \note This constructor is not implemented and should not be used.
   */
  FeedbackForcer();

  /*!
   * \brief Copy constructor.
   *
   * \note This constructor is not implemented and should not be used.
   *
   * \param from The value to copy to this object.
   */
  FeedbackForcer(const FeedbackForcer& from);

  /*!
   * \brief Assignment operator.
   *
   * \note This operator is not implemented and should not be used.
   *
   * \param that The value to assign to this object.
   *
   * \return A reference to this object.
   */
  FeedbackForcer& operator=(const FeedbackForcer& that);

  const VelocityBcCoefs* d_velocity_bc;

  const INSHierarchyIntegrator* const d_fluid_solver;
  Pointer<PatchHierarchy<NDIM> > d_patch_hierarchy;
  
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <FeedbackForcer.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FeedbackForcer
