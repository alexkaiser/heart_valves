// Filename: FeedbackForcer.cpp
// Created on 04 May 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

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



// #define FLOW_STRAIGHTENER
#define OPEN_BOUNDARY_STABILIZATION

#define FLOW_AVERAGER

#define FULL_FLOW_CLAMP 
#define FULL_FLOW_CLAMP_TIME 0.01

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
                               CirculationModel_with_lv* circ_model_with_lv,
                               CirculationModel_RV_PA* circ_model_rv_pa) 
  : d_fluid_solver(fluid_solver), 
    d_patch_hierarchy(patch_hierarchy),
    d_circ_model_with_lv(circ_model_with_lv), 
    d_circ_model_rv_pa(circ_model_rv_pa)
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
    int axis = 2;
    
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
        if (d_circ_model_with_lv){

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
                
                // always working on the z component, no need to 
                int component = 2;
                
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++){

                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                    
                    double X[NDIM];

                    for (int d = 0; d < NDIM; ++d){
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                    }

                    const int in_aorta  = d_circ_model_with_lv->point_in_aorta(X[0],X[1]); 
                    const int in_atrium = d_circ_model_with_lv->point_in_atrium(X[0],X[1]); 

                    if (in_aorta && in_atrium){
                        TBOX_ERROR("Position is within both aorta and atrium, should be impossible\n"); 
                    }

                    // no bdry stab unless one of the conditionals is met 
                    double mask = 0.0;

                    double U_goal = 0.0; 

                    if (in_aorta){
                        const double n = is_lower ? -1.0 : +1.0;
                        const double U_dot_n = U * n;
                        if ((inflow_bdry_aorta && U_dot_n > 0.0) || (outflow_bdry_aorta && U_dot_n < 0.0)){
                            mask = 1.0;
                        }

                        #ifdef FLOW_AVERAGER
                            // set goal to be equal to average flow 
                            if (d_circ_model_with_lv->d_area_initialized){
                                U_goal = d_circ_model_with_lv->d_Q_aorta / d_circ_model_with_lv->d_area_aorta;
                                mask = 1.0;
                                //pout << "In averager woot, aorta flux = " <<  d_circ_model_with_lv->d_Q_aorta << ", aorta area = " << d_circ_model_with_lv->d_area_aorta << "mean flow aorta = " << U_goal << "\n"; 
                            }
                        #endif

                    }

                    if (in_atrium){
                        const double n = is_lower ? -1.0 : +1.0;
                        const double U_dot_n = U * n;
                        if ((inflow_bdry_left_atrium && U_dot_n > 0.0) || (outflow_bdry_left_atrium && U_dot_n < 0.0)){
                            mask = 1.0;
                        }

                        #ifdef FLOW_AVERAGER
                            // set goal to be equal to average flow 
                            if (d_circ_model_with_lv->d_area_initialized){
                                U_goal = d_circ_model_with_lv->d_Q_left_atrium / d_circ_model_with_lv->d_area_atrium; 
                                mask = 1.0;
                                //pout << "In averager woot, atrium flux = " <<  d_circ_model_with_lv->d_Q_left_atrium << ", atrium area = " << d_circ_model_with_lv->d_area_atrium << "mean flow atrium = " << U_goal << "\n"; 
                            }
                        #endif

                    }

                    if (mask > 0.0){
                        const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                        mask *= smooth_kernel((X[axis] - x_bdry) / L);
                        (*F_data)(i_s) += mask * (-kappa * (U - U_goal));
                    }
                }
                
            }
        } // if (d_circ_model_with_lv)


        if (d_circ_model_rv_pa){

            // Attempt to prevent flow reversal points near the domain boundary.
            double height_physical = 0.2;

            // hardcoded z axis top here 
            for(int axis=0; axis<NDIM; axis++){ 
                for(int side=0; side<2; side++){ 

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
                        
                        // here check all components 
                        for (int component = 0; component<NDIM; component++){
                        
                            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++){

                                const Index<NDIM>& i = b();
                                const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                                const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                                const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                                const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                                
                                double X[NDIM];

                                for (int d = 0; d < NDIM; ++d){
                                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                                }

                                double X_in_plane_1 = 0.0; 
                                double X_in_plane_2 = 0.0; 
                                if (axis == 0)
                                {
                                    X_in_plane_1 = X[1]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 1)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 2)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[1]; 
                                }
                                else{
                                    TBOX_ERROR("Invalid axis\n"); 
                                }

                                const int in_right_ventricle  = d_circ_model_rv_pa->point_in_right_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_right_pa         = d_circ_model_rv_pa->point_in_right_pa       (X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_left_pa          = d_circ_model_rv_pa->point_in_left_pa        (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_right_ventricle && in_right_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }
                                if (in_right_ventricle && in_left_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }
                                if (in_right_pa && in_left_pa){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }

                                // no bdry stab unless one of the conditionals is met 
                                double mask = 0.0;
                                double U_goal = 0.0; 

                                if (in_right_ventricle){

                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_right_ventricle / d_circ_model_rv_pa->d_area_right_ventricle;
                                            }
                                            mask = 1.0;
                                        }
                                    #endif

                                }

                                if (in_right_pa){


                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_right_pa / d_circ_model_rv_pa->d_area_right_pa;
                                            }
                                            mask = 1.0;
                                        }
                                    #endif

                                }

                                if (in_left_pa){

                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_left_pa / d_circ_model_rv_pa->d_area_left_pa;
                                            }
                                            mask = 1.0;
                                        }
                                    #endif

                                }

                                if (mask > 0.0){
                                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                                    mask *= smooth_kernel((X[axis] - x_bdry) / (height_physical/2.0));
                                    (*F_data)(i_s) += mask * (-kappa * (U - U_goal));
                                }

                            }
                        }
                        
                    }
                } // for(int side=0; side<2; side++)
            } // for(int axis=0; axis<NDIM; axis++)
        } // if (d_circ_model_rv_pa)


    #endif 
    
    #ifdef FULL_FLOW_CLAMP
        if (data_time < FULL_FLOW_CLAMP_TIME){
            F_data->fillAll(0.0);

            // linear decrease in coefficient value 
            // from max 
            double k_full_clamp; 
            if (cycle_num > 0){
                k_full_clamp = (1 - data_time/FULL_FLOW_CLAMP_TIME) * 0.25 * rho / dt;
            }
            else{
                k_full_clamp = 0.0; 
            }

            // Clamp the velocity in all components
            for (int component = 0; component < NDIM; ++component){
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){

                    const Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);

                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new     = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U         = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;

                    (*F_data)(i_s)        += -k_full_clamp * U;
                }
            }
        }
    #endif


    return;
} // setDataOnPatch


/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ////////////
