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
#define OPEN_BOUNDARY_STABILIZATION

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

FeedbackForcer::FeedbackForcer(const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               CirculationModel_with_lv* circ_model_with_lv) 
  : d_fluid_solver(fluid_solver), 
    d_patch_hierarchy(patch_hierarchy),
    d_circ_model_with_lv(circ_model_with_lv)
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
    
//        const double* const x_lower_global = grid_geometry->getXLower();
        const double* const upper_limit_global = grid_geometry->getXUpper();

        // physical height of region of stabilization
        // stabilization is smoothed out from bottom of domain to here
        const double height_physical = 1.0;
        const double min_height_force_applied = upper_limit_global[axis] - height_physical;        
        
        double k_straightener[NDIM];
        k_straightener[0] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
        k_straightener[1] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
        k_straightener[2] = cycle_num >= 0 ?             0.0 : 0.0; // no friction in the z direction
                                                                   
        // Clamp the velocity in the x,y components
        for (int component = 0; component < NDIM; ++component){
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){
            
                const Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);

                // get the height, which determines whether there is force
                const double z = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                
                if (z > min_height_force_applied){
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new     = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U         = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;

                    const double weight    = smooth_kernel((z - upper_limit_global[axis]) / height_physical);
                
                    (*F_data)(i_s)        += weight*(-k_straightener[component] * U);
                }
            }
        }
    
    #endif
    
    #ifdef OPEN_BOUNDARY_STABILIZATION

        // Attempt to prevent flow reversal points near the domain boundary.

        // hardcoded z axis top here 
        // static const int axis = 2;
        int side = 1; 

        const double L = max(dx_coarsest[axis], 2.0 * dx_finest[axis]);
        const int offset = static_cast<int>(L / dx[axis]);
        
        const bool is_lower = side == 0;

        const bool inflow_bdry_aorta        = (d_circ_model_with_lv->d_Q_aorta < 0.0);
        const bool outflow_bdry_aorta       = !(inflow_bdry_aorta);
        const bool inflow_bdry_left_atrium  = (d_circ_model_with_lv->d_Q_left_atrium < 0.0);
        const bool outflow_bdry_left_atrium = !(inflow_bdry_left_atrium);

        if (pgeom->getTouchesRegularBoundary(axis, side)){
            Box<NDIM> bdry_box = domain_box;
            if (side == 0){
                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
            }
            else{
                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
            }
            bdry_box = bdry_box * patch_box;
            for (int component = 0; component < NDIM; ++component){
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++){

                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    
                    double X[NDIM];
                    double dist_sq_aorta = 0.0;
                    double dist_sq_atrium = 0.0;

                    for (int d = 0; d < NDIM; ++d){
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                        if (d != axis){
                            dist_sq_aorta  += pow(X[d] - d_circ_model_with_lv->d_center_aorta[d],  2.0);
                            dist_sq_atrium += pow(X[d] - d_circ_model_with_lv->d_center_atrium[d], 2.0);
                        }
                    }

                    const double dist_aorta  = sqrt(dist_sq_aorta);
                    const double dist_atrium = sqrt(dist_sq_atrium);

                    double mask = (component == axis ? 0.0 : 1.0);

                    if (component == axis && (dist_aorta < d_circ_model_with_lv->d_radius_aorta)){
                        const double n = is_lower ? -1.0 : +1.0;
                        const double U_dot_n = U * n;
                        if ((inflow_bdry_aorta && U_dot_n > 0.0) || (outflow_bdry_aorta && U_dot_n < 0.0)){
                            mask = 1.0;
                        }
                    }

                    if (component == axis && (dist_atrium < d_circ_model_with_lv->d_radius_atrium)){
                        const double n = is_lower ? -1.0 : +1.0;
                        const double U_dot_n = U * n;
                        if ((inflow_bdry_left_atrium && U_dot_n > 0.0) || (outflow_bdry_left_atrium && U_dot_n < 0.0)){
                            mask = 1.0;
                        }
                    }

                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    mask *= smooth_kernel((X[axis] - x_bdry) / L);
                    (*F_data)(i_s) += mask * (-kappa * U);
                }
            }
        }
        


        // Open boundary stabilization
/*        double width[NDIM];
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
                    
                    // If hitting instabilities, clobber all other forcing values 
                    if ((inflow_bdry && (U * n > 0.0)) || (outflow_bdry && (U * n < 0.0))){
                        const double x = x_lower[axis] + dx[axis] * static_cast<double>(i(axis) - patch_box.lower(axis));
                        const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                        (*F_data)(i_s) = smooth_kernel((x - x_bdry) / width[axis]) * kappa * (0.0 - U);
                    }
                }
            } // if outflow_bdry
        } // side */ 
    #endif 
    
    return;
} // setDataOnPatch


/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ////////////
