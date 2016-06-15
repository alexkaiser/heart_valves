
// FourierBodyForce.cpp
// Created Alex Kaiser, 6/2016


// Modified from:
// Filename: FeedbackForcer.cpp
// Created on 04 May 2007 by Boyce Griffith

#include "FourierBodyForce.h"

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
#include <SideData.h>
#include <tbox/Utilities.h>



/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/*
namespace
{
inline double
smooth_kernel(const double r)
{
  return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
}
*/


////////////////////////////// PUBLIC ///////////////////////////////////////

FourierBodyForce::FourierBodyForce(const fourier_series_data* fourier,
                               const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
  : d_fourier(fourier), d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
{
  // intentionally blank
  return;
} // FourierBodyForce

FourierBodyForce::~FourierBodyForce()
{
  // intentionally blank
  return;
} // ~FourierBodyForce

bool
FourierBodyForce::isTimeDependent() const
{
  return true;
} // isTimeDependent

void
FourierBodyForce::setDataOnPatch(const int data_idx,
                               Pointer<Variable<NDIM> > /*var*/,
                               Pointer<Patch<NDIM> > patch,
                               const double data_time,
                               const bool initial_time,
                               Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
    #if !defined(NDEBUG)
        TBOX_ASSERT(F_data);
    #endif
    
    F_data->fillAll(0.0);
    
    //if (initial_time) return;
    
    // get basic information
//    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
//    const double dt = d_fluid_solver->getCurrentTimeStepSize();
//    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    
    // this is a constant but should be messed with
    //const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    
    // pull velocity variables
/*    Pointer<SideData<NDIM, double> > U_current_data =
    patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    Pointer<SideData<NDIM, double> > U_new_data =
    patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());

    #if !defined(NDEBUG)
        TBOX_ASSERT(U_current_data);
    #endif
*/

    
    // stuff about the physical box and mesh structures
    // take as needed
    const Box<NDIM>& patch_box = patch->getBox();
    
/*    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    */
    
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    
    // Added these
    const double* const x_lower = grid_geometry->getXLower();
    const double* const x_upper = grid_geometry->getXUpper();
    std::cout << "code thinks the lower domain is " << x_lower[0] << ", " << x_lower[0] << ", " << x_lower[0] << "\n" ;
    std::cout << "code thinks the upper domain is " << x_upper[0] << ", " << x_upper[0] << ", " << x_upper[0] << "\n" ;
    
    const double z_domain_length = x_upper[2] - x_lower[2];
    
    /*
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());
*/
    
    // reduce the current beat modulo the fourier bc
    //const double beat_time = (d_velocity_bc->d_fourier)->L;
    //const double t_reduced = data_time - floor(data_time / beat_time);
    
    // index without periodicity
    unsigned int k = (unsigned int) floor(data_time / (d_fourier->dt));
    
    // take periodic reduction                         
    unsigned int idx = k % (d_fourier->N_times);
    
    
    const double force = -MMHG_TO_CGS * d_fourier->values[idx] / z_domain_length;

    // Always force in the negative z direction
    const int component = 2;
    //for (int component = 0; component < NDIM; ++component)
    
    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){
        const Index<NDIM>& i = b();
        const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);

        (*F_data)(i_s) = force; // FORCE GOES HERE
    }
    
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ////////////
