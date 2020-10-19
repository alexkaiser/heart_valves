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

#include <queue>

// #define FLOW_STRAIGHTENER
#define OPEN_BOUNDARY_STABILIZATION

// #define FLOW_AVERAGER

#define FULL_FLOW_CLAMP 
#define FULL_FLOW_CLAMP_TIME 0.1

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
                               CirculationModel_with_lv* circ_model_with_lv, // = NULL 
                               CirculationModel_RV_PA* circ_model_rv_pa, // = NULL 
                               CirculationModel_aorta* circ_model_aorta, // = NULL 
                               bool damping_outside, // = false 
                               string lag_file_name, // = ""
                               string internal_ring_file_name) // = "" 
  : d_fluid_solver(fluid_solver), 
    d_patch_hierarchy(patch_hierarchy),
    d_circ_model_with_lv(circ_model_with_lv), 
    d_circ_model_rv_pa(circ_model_rv_pa),
    d_circ_model_aorta(circ_model_aorta),
    d_damping_outside(damping_outside), 
    d_damping_initialized(false)
{

    if(d_damping_outside){

        if (lag_file_name.empty()){
            TBOX_ERROR("Must provide valid lag_file_name if damping_outside is on");
        }
        if (internal_ring_file_name.empty()){
            TBOX_ERROR("Must provide valid internal_ring_file_name if damping_outside is on");
        }

        this->initialize_masks(d_patch_hierarchy->getGridGeometry(),
                                              lag_file_name,
                                              internal_ring_file_name);
    }


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
    

    if (d_damping_outside){
        
        if(!d_damping_initialized){
            TBOX_ERROR("damping on but not initialized"); 
        }

        // Clamp the velocity in the interior of the rigid solid.
        for (int component = 0; component < NDIM; ++component){
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){
                const Index<NDIM>& i = b();
                const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;
                
                double X[NDIM];
                for (int d = 0; d < NDIM; ++d){
                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                }

                int idx = this->get_1d_idx(X); 

                // top mask may be 1/2 mesh width over because of staggered scheme 
                // just skip these points 
                double mask = (idx < d_N_Eulerian_total) ? d_masks_linear_array[idx] : 0.0; 

                (*F_data)(i_s) += (mask) * (-kappa * U);

                // if ((idx < 0) || (idx >= d_N_Eulerian_total)){
                //     pout << "bad index, idx = " << idx << ", X = " << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
                //     pout << "in outside damping, (*F_data)(i_s) = " << (*F_data)(i_s) << ", U = " <<  U << ", mask = " << mask << "\n\n"; 
                // }
                if (!((mask == 0.0) || (mask == 1.0))){
                    pout << "bad mask value without bad index:\n"; 
                    pout << "in outside damping, (*F_data)(i_s) = " << (*F_data)(i_s) << ", U = " <<  U << ", mask = " << mask << "\n\n"; 
                }

            }
        }
    }



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
                    double mask = 1.0;

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

                                // bool tangential_damp_to_zero = false;

                                if (in_right_ventricle){

                                    const double n = is_lower ? -1.0 : +1.0;
                                    const double U_dot_n = U * n;

                                    if (axis == component){
                                        // different signs give local flow reversal
                                        if ((d_circ_model_rv_pa->d_Q_right_ventricle * U_dot_n) < 0.0){
                                            mask = 1.0;
                                        }
                                    }

                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_right_ventricle / d_circ_model_rv_pa->d_area_right_ventricle;
                                                mask = 1.0;
                                            }
                                        }
                                    #endif

                                    // if ((axis != component) && (tangential_damp_to_zero)){
                                    //     mask = 1.0;
                                    // }

                                }

                                if (in_right_pa){

                                    const double n = is_lower ? -1.0 : +1.0;
                                    const double U_dot_n = U * n;

                                    if (axis == component){
                                        // different signs give local flow reversal
                                        if ((d_circ_model_rv_pa->d_Q_right_pa * U_dot_n) < 0.0){
                                            mask = 1.0;
                                        }
                                    }

                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_right_pa / d_circ_model_rv_pa->d_area_right_pa;
                                                mask = 1.0;
                                            }
                                        }
                                    #endif

                                    // if ((axis != component) && (tangential_damp_to_zero)){
                                    //     mask = 1.0;
                                    // }

                                }

                                if (in_left_pa){

                                    const double n = is_lower ? -1.0 : +1.0;
                                    const double U_dot_n = U * n;

                                    if (axis == component){
                                        // different signs give local flow reversal
                                        if ((d_circ_model_rv_pa->d_Q_left_pa * U_dot_n) < 0.0){
                                            mask = 1.0;
                                        }
                                    }

                                    #ifdef FLOW_AVERAGER
                                        // set goal to be equal to average flow 
                                        if (d_circ_model_rv_pa->d_area_initialized){
                                            if (axis == component){
                                                U_goal = d_circ_model_rv_pa->d_Q_left_pa / d_circ_model_rv_pa->d_area_left_pa;
                                                mask = 1.0;
                                            }
                                        }
                                    #endif

                                    // if ((axis != component) && (tangential_damp_to_zero)){
                                    //     mask = 1.0;
                                    // }

                                }

                                if (mask > 0.0){
                                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                                    mask *= smooth_kernel((X[axis] - x_bdry) / (height_physical/2.0));
                                    (*F_data)(i_s) += mask * (-kappa * (U - U_goal));
                                }

                            }
                        } // for (int component = 0; component<NDIM; component++)
                    } // if (pgeom->getTouchesRegularBoundary(axis, side)){
                } // for(int side=0; side<2; side++)
            } // for(int axis=0; axis<NDIM; axis++)
        } // if (d_circ_model_rv_pa)



        if (d_circ_model_aorta){

            // Attempt to prevent flow reversal points near the domain boundary.
            double height_physical = 0.2;

            // hardcoded z axis top here 
            for(int axis=0; axis<NDIM; axis++){ 
                for(int side=0; side<2; side++){ 

                    const double L = max(dx_coarsest[axis], 2.0 * dx_finest[axis]);
                    const int offset = static_cast<int>(L / dx[axis]);
                    const bool is_lower = side == 0;

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

                                const int in_ventricle  = d_circ_model_aorta->point_in_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_aorta      = d_circ_model_aorta->point_in_aorta    (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_ventricle && in_aorta){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }

                                // no bdry stab unless one of the conditionals is met 
                                double mask = 0.0;
                                double U_goal = 0.0; 

                                bool tangential_damp_to_zero = true;

                                if (in_ventricle){

                                    const double n = is_lower ? -1.0 : +1.0;
                                    const double U_dot_n = U * n;

                                    if (axis == component){
                                        // different signs give local flow reversal
                                        if ((d_circ_model_aorta->d_Q_ventricle * U_dot_n) < 0.0){
                                            mask = 1.0;
                                        }
                                    }

                                    // #ifdef FLOW_AVERAGER
                                    // set goal to be equal to average flow 
                                    if (d_circ_model_aorta->d_area_initialized){
                                        if (axis == component){
                                            U_goal = d_circ_model_aorta->d_Q_ventricle / d_circ_model_aorta->d_area_ventricle;
                                            mask = 1.0;
                                        }
                                    }
                                    // #endif

                                    if ((axis != component) && (tangential_damp_to_zero)){
                                        mask = 1.0;
                                    }

                                }
                                else if (in_aorta){

                                    const double n = is_lower ? -1.0 : +1.0;
                                    const double U_dot_n = U * n;

                                    if (axis == component){
                                        // different signs give local flow reversal
                                        if ((d_circ_model_aorta->d_Q_aorta * U_dot_n) < 0.0){
                                            mask = 1.0;
                                        }
                                    }

                                    // #ifdef FLOW_AVERAGER
                                    //     // set goal to be equal to average flow 
                                    //     if (d_circ_model_aorta->d_area_initialized){
                                    //         if (axis == component){
                                    //             U_goal = d_circ_model_aorta->d_Q_aorta / d_circ_model_aorta->d_area_aorta;
                                    //             mask = 1.0;
                                    //         }
                                    //     }
                                    // #endif

                                    // if ((axis != component) && (tangential_damp_to_zero)){
                                    //     mask = 1.0;
                                    // }

                                }
                                else{
                                    // damp all edges outside inlets and outlets
                                    mask = 1.0; 
                                }

                                if (mask > 0.0){
                                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                                    mask *= smooth_kernel((X[axis] - x_bdry) / (height_physical/2.0));
                                    (*F_data)(i_s) += mask * (-kappa * (U - U_goal));
                                }

                            }
                        } // for (int component = 0; component<NDIM; component++)
                    } // if (pgeom->getTouchesRegularBoundary(axis, side)){
                } // for(int side=0; side<2; side++)
            } // for(int axis=0; axis<NDIM; axis++)
        } // if (d_circ_model_aorta)

    #endif 
    
    #ifdef FULL_FLOW_CLAMP
        if (data_time < FULL_FLOW_CLAMP_TIME){
            // F_data->fillAll(0.0);

            // linear decrease in coefficient value 
            // from max 
            double k_full_clamp; 
            if (cycle_num > 0){
                // k_full_clamp = (1 - data_time/FULL_FLOW_CLAMP_TIME) * 0.25 * rho / dt;
                k_full_clamp = 0.25 * rho / dt;
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

                    if ((*F_data)(i_s) == 0.0){
                        (*F_data)(i_s) += -k_full_clamp * U;
                    }

                }
            }
        }
    #endif


    return;
} // setDataOnPatch

void FeedbackForcer::initialize_masks(Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
                                      string lag_file_name,
                                      string internal_ring_file_name){
    
    // Serial code, redundant in parallel 
    //
    // Marks all points outside of current Lagrangian stucture 
    // 
        
    // set up the data
    
    // pout << "got some data from the ldata manager\n";
    
    // get data from petsc arrays
    int N_Lagrangian;            // keeps the full list of all three components 

    // read file here 
    ifstream f; 
    f.open(lag_file_name.c_str(), ios::in);
    
    if (!f.is_open()){
        std::cout << "Failed to open file\n" ; 
    } 
    
    // total lagrangian points 
    f >> N_Lagrangian;

    // flattened array of lagrangian points 
    double *X_lagrangian = new double[3*N_Lagrangian]; 
    for (int j=0; j<(3*N_Lagrangian) && (!f.eof()); j++){
        f >> X_lagrangian[j]; 
    }

    f.close(); 

    // read file for seed point 
    ifstream f_internal_ring; 
    f_internal_ring.open(internal_ring_file_name.c_str(), ios::in);
    
    if (!f_internal_ring.is_open()){
        std::cout << "Failed to open file\n" ; 
    } 
    
    // lagrangian 
    int N_ring; 
    double internal_point[3] = {0.0, 0.0, 0.0}; 

    double x_tmp, y_tmp, z_tmp; 

    f_internal_ring >> N_ring;

    for (int j=0; j<N_ring && (!f_internal_ring.eof()); j++){
        f_internal_ring >> x_tmp; 
        f_internal_ring >> y_tmp; 
        f_internal_ring >> z_tmp; 
        internal_point[0] += x_tmp; 
        internal_point[1] += y_tmp; 
        internal_point[2] += z_tmp;         
    }
    internal_point[0] /= N_ring; 
    internal_point[1] /= N_ring; 
    internal_point[2] /= N_ring; 

    f_internal_ring.close();     

    // set geometry information 
    const double *dx_each = grid_geometry->getDx(); 
    d_bdry_low            = grid_geometry->getXLower(); 
    d_bdry_up             = grid_geometry->getXUpper(); 
    
    // consistency check 
    if ( (dx_each[0] - dx_each[1] > 1e-12) || (dx_each[0] - dx_each[2] > 1e-12) ){
        pout << "All meshes must be the same size to compute Lagrangian error\n" ; 
        pout << "Mesh dimensions = (" << dx_each[0] << ", " << dx_each[1] << ", " << dx_each[2] << ")\n"; 
        abort(); 
    } 

    d_dx = dx_each[0]; 
        
    // unsigned int N[NDIM]; 
    for(int i=0; i<NDIM; i++){
        d_N[i] = (int) ((d_bdry_up[i] - d_bdry_low[i]) / d_dx); 
    }
    
    d_N_Eulerian_total = d_N[0] * d_N[1] * d_N[2]; 

    int *indices_one_dimensional = new int[d_N_Eulerian_total]; 
    d_masks_linear_array = new double[d_N_Eulerian_total]; 
    
    for(int i=0; i<d_N_Eulerian_total; i++){
        indices_one_dimensional[i] = 0; 
    }

    int idx[3]; 
    int idx_nbr[3];             // may be negative but then this is out of bounds 
    int idx_one_dimensional; 
    int idx_nbr_one_dimensional; 
        
    // Mark all elements near Lagrangian boundary. 
    // This marks 8 points on the edges of the cube which contains the Lagrangian point. 
    // Lower boundary is inclusive, upper exclusive. 
    //     
    int eulerian_pts_marked = 0; 
    
    for(int idx_lag=0; idx_lag<(3*N_Lagrangian); idx_lag+=3){
        
        idx_one_dimensional = this->get_1d_idx(X_lagrangian+idx_lag); 
        
        this->get_3d_idx(idx_one_dimensional, idx); 
        
        int min_range = -2; 
        int max_range = 2; 

        for(int i=min_range; i<=max_range; i++){
            for(int j=min_range; j<=max_range; j++){
                for(int k=min_range; k<=max_range; k++){
                
                    idx_nbr[0] = idx[0] + i; 
                    idx_nbr[1] = idx[1] + j; 
                    idx_nbr[2] = idx[2] + k;                     

                    if( this->in_bounds(idx_nbr) ){
                        
                        idx_nbr_one_dimensional = this->get_1d_idx_from_3d_idx(idx_nbr); 
                        
                        // mark it 
                        if (indices_one_dimensional[idx_nbr_one_dimensional] == 0){
                            indices_one_dimensional[idx_nbr_one_dimensional] = 2;
                            eulerian_pts_marked++;                             
                        }
                         
                    }  
                }
            }
        }
    }


    std::queue<int> indices_queue; 
       
    // physical points to 1d index 
    idx_one_dimensional = this->get_1d_idx(internal_point); 
   
    pout << "idx_one_dimensional = " << idx_one_dimensional << "\n"; 
    pout << "internal_point = " << internal_point[0] << ", " << internal_point[1] << ", " << internal_point[2] << "\n"; 

    if ((idx_one_dimensional < 0) || (idx_one_dimensional > d_N_Eulerian_total)){
        TBOX_ERROR("initial index out of range");
    }

    if(indices_one_dimensional[idx_one_dimensional] != 0){
        TBOX_ERROR("index of internal point is aleardy marked \n");
    }
    
    indices_queue.push( idx_one_dimensional ); 
    indices_one_dimensional[idx_one_dimensional] = 1; 
    
    int total_internal = 0; 
 
    while( !indices_queue.empty() ){
        
        // got a point... 
        total_internal++; 
        
        idx_one_dimensional = indices_queue.front();
        indices_queue.pop();

        this->get_3d_idx(idx_one_dimensional, idx); 

        for(int i=-1; i<=1; i++){
            for(int j=-1; j<=1; j++){
                for(int k=-1; k<=1; k++){
                                        
                    // don't compare with self
                    if((i==0) && (j==0) && (k==0))
                        continue; 
                
                    idx_nbr[0] = idx[0] + i; 
                    idx_nbr[1] = idx[1] + j; 
                    idx_nbr[2] = idx[2] + k;                     
                                    
                    // is this a vaild index? 
                    if( this->in_bounds(idx_nbr) ){
                                                
                        idx_nbr_one_dimensional = this->get_1d_idx_from_3d_idx(idx_nbr); 
                        
                        // is it unfilled? 
                        if(indices_one_dimensional[idx_nbr_one_dimensional] == 0){
                        
                            // mark it 
                            indices_one_dimensional[idx_nbr_one_dimensional] = 1; 
                            
                            // push on 
                            indices_queue.push(idx_nbr_one_dimensional); 
                        }
                    }  
                }
            }
        }
    }
 
    // finally set the masks 
    for (int j=0; j<d_N_Eulerian_total; j++){
        if (indices_one_dimensional[j] == 0){
            // full mask 
            d_masks_linear_array[j] = 1.0; 
        }
        else{
            // no mask, on Lag structure 
            // or inside of it 
            d_masks_linear_array[j] = 0.0; 
        }

    }

    bool debug_plot_file = false; 
    if (debug_plot_file){
        std::ofstream mask_data;
        mask_data.open("mask_data.csv", ios_base::out | ios_base::trunc);

        mask_data << "x, y, z, v \n"; 

        for(int i=0; i<d_N_Eulerian_total; i++){
            get_3d_idx(i, idx); 

            double x = d_dx * idx[0] + d_bdry_low[0]; 
            double y = d_dx * idx[1] + d_bdry_low[1]; 
            double z = d_dx * idx[2] + d_bdry_low[2]; 

            if (indices_one_dimensional[i] != 0){
                mask_data << x << ", " << y << ", " << z << ", " << indices_one_dimensional[i] << "\n";                 
            }


        }
        mask_data.close(); 

    }


    delete[] indices_one_dimensional; 

    bool debug_vol = true;  
    if (debug_vol){
        double internal_vol = ((double) total_internal) * d_dx*d_dx*d_dx;
        double total_vol = ((double) total_internal + eulerian_pts_marked) * d_dx*d_dx*d_dx;
        double total_eulerian_vol = d_dx*d_dx*d_dx * d_N[0] * d_N[1] * d_N[2];
        pout << "found intenal volume " << internal_vol << ", internal and lag volume " << total_vol << ", total eulerian vol = " << total_eulerian_vol << "\n";
    }

    d_damping_initialized = true; 

}


inline unsigned int FeedbackForcer::get_1d_idx(const double *point){
    /*
    Computes the global 1-d index of the given point. 
    Assumes x,y,z order.  
    
    Input: 
    const double *point      Cartesian coordinates of the point. 
    const double *bdry_low   Lower boundary of domain 
    const int *N             Number points in domain 
    const double dx          Mesh width
    
    Output 
    int (returned)           Global index 
    */

    int idx_x = (int) ((point[0] - d_bdry_low[0]) / d_dx); 
    int idx_y = (int) ((point[1] - d_bdry_low[1]) / d_dx);         
    int idx_z = (int) ((point[2] - d_bdry_low[2]) / d_dx); 
    
    return idx_x + idx_y * d_N[0] + idx_z * d_N[0] * d_N[1];
}


inline unsigned int FeedbackForcer::get_1d_idx_from_3d_idx(const int *idx){
    /*
    Computes the global 1-d index from the 3 variable local index
    Assumes x,y,z order.  
    
    Input: 
    const int *N               Number points in domain 
    const unsigned int *idx    Index in 3 coordinates  
    
    Output 
    int (returned)             Global index 
    */

    int idx_x = (unsigned int) idx[0]; 
    int idx_y = (unsigned int) idx[1]; 
    int idx_z = (unsigned int) idx[2]; 

    return idx_x + idx_y * d_N[0] + idx_z * d_N[0] * d_N[1];
}


inline void FeedbackForcer::get_3d_idx(const int one_dimensional_idx, int *idx){
    /*
    Computes the 3 coordinates from the one idex
    Assumes x,y,z order.  
    
    Input: 
    const unsigned int one_dimensional_idx      Cartesian coordinates of the point. 
    const int *N                       Number points in domain 

    Output:
    unsigned int *idx                  Three coordinates of the index. 
    */

    unsigned int temp; 
    unsigned int i,j,k;
    unsigned int one_dimensional_idx_unsigned = (unsigned int) one_dimensional_idx; 

    i      = one_dimensional_idx_unsigned % d_N[0];           // first index from modding out the first dimension 
    temp   = (one_dimensional_idx_unsigned-i) / d_N[0];       // subtract from global and divide out first, gives   temp = j + k*N[1]  
    j      = temp % d_N[1];                 // get j 
    k      = (temp - j)/d_N[1];             // temp - j = k*N[1]
    
    idx[0] = i; 
    idx[1] = j; 
    idx[2] = k;     
}


inline bool FeedbackForcer::in_bounds(const int *idx){
    /*
    Checks whether the index is in bounds or not. 
    
    Input: 
    const int *N                       Number points in domain 
    unsigned int *idx                  Three coordinates of the index. 
    
    Output:
    True if in bounds. 
    */

    if((idx[0] < 0) || (idx[1] < 0) || (idx[2] < 0))
        return false; 
        
    if((idx[0] >= ((int) d_N[0]))  ||  (idx[1] >= ((int) d_N[1]))  ||  (idx[2] >= ((int) d_N[2]))) 
        return false; 
        
    return true; 
    
}


/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ////////////
