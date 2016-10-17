// Filename: FeedbackForcer.cpp
// Created on 04 May 2007 by Boyce Griffith

#include "FeedbackForcer.h"

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



#define FLOW_STRAIGHTENER

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline double
smooth_kernel(const double r)
{
  return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
}

////////////////////////////// PUBLIC ///////////////////////////////////////

FeedbackForcer::FeedbackForcer(const VelocityBcCoefs* velocity_bc,
                               const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
  : d_velocity_bc(velocity_bc), d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy)
{
  // intentionally blank
  return;
} // FeedbackForcer

FeedbackForcer::~FeedbackForcer()
{
  // intentionally blank
  return;
} // ~FeedbackForcer

bool
FeedbackForcer::isTimeDependent() const
{
  return true;
} // isTimeDependent

void
FeedbackForcer::setDataOnPatch(const int data_idx,
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
    
    if (initial_time) return;
    
    // get basic information
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    
    // this is a constant but should be messed with
    const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    
    // pull velocity variables
    Pointer<SideData<NDIM, double> > U_current_data =
    patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    Pointer<SideData<NDIM, double> > U_new_data =
    patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(U_current_data);
    #endif

    
    // stuff about the physical box and mesh structures
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());
    
    
    // hardcoded for z axis for now
    static const int axis = 2;
    
    // Flow straightener and friction (if desired)
    #ifdef FLOW_STRAIGHTENER
    
        const double* const x_lower_global = grid_geometry->getXLower();
        
        // physical height of region of stabilization
        // stabilization is smoothed out from bottom of domain to here
        const double height_physical = 1.0;
        const double max_height_force_applied = height_physical + x_lower_global[2];
        const double center = x_lower_global[axis] + 0.5*height_physical;
        
        // this may be very, very large
        // consider changing it
        double k_straightener[NDIM];
        k_straightener[0] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
        k_straightener[1] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
        k_straightener[2] = cycle_num >= 0 ?           0.0e4 : 0.0; // much lower friction in the z direction
                                                           // at U = 10cm/s, this is ~10x force of gravity 
        
        // Clamp the velocity in the x,y components
        // Clamp the velocity in the z component, but a lot less
        for (int component = 0; component < NDIM; ++component){
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){
            
                const Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);

                // get the height, which determines whether there is force
                const double z = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                
                if (z < max_height_force_applied){
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new     = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U         = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;

                    const double weight    = smooth_kernel((z - center) / (dx[axis]*height_physical));
                
                    (*F_data)(i_s)        += weight*(-k_straightener[component] * U);
                }
            }
        }
    
    #endif
    
    
    // Open boundary stabilization
    double width[NDIM];
    for(int i=0; i<NDIM; i++){
        width[i] = 4.0 * dx_finest[i];
    }

    // Attempt to prevent flow reversal points near the domain boundary.
    for (int side = 0; side <= 1; ++side){
        const bool is_lower = (side == 0);
        
        // Upper z face outward flux
        const double qsrc = (d_velocity_bc->d_circ_model)->d_qsrc[0];
        
        // Upper z direction in and out flows
        bool inflow_bdry  = (qsrc < 0.0);
        bool outflow_bdry = (qsrc > 0.0);

        // Reverse signs for lower boundaries
        if (is_lower){
            inflow_bdry  = (!inflow_bdry);
            outflow_bdry = (!outflow_bdry);
        }
        
        // don't need to mess with this loop if it is not an outflow boundary
        
        // pulled directly from boundary stab code
        // only change is to add "outflow boundary" at all relevant places
        
        if (pgeom->getTouchesRegularBoundary(axis, side)){
            
            Box<NDIM> bdry_box = domain_box;
            
            // changed this to use width (rather than boundary stabilizer instance var d_width)
            const int offset = static_cast<int>(width[axis] / dx[axis]);
            
            if (is_lower){
                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
            }
            else{
                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
            }
            
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box * patch_box, axis)); b; b++){
                const Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                const double n = is_lower ? -1.0 : +1.0;
                
                if ((inflow_bdry && (U * n > 0.0)) || (outflow_bdry && (U * n < 0.0))){
                    const double x = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    (*F_data)(i_s) = smooth_kernel((x - x_bdry) / width[axis]) * kappa * (0.0 - U);
                }
            }
        } // if outflow_bdry
    } // side
    
    return;
} // setDataOnPatch


/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ////////////
