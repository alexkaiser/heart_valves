
// FourierBodyForce.cpp
// Created Alexander D. Kaiser, 6/2016


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


#define FLOW_STRAIGHTENER

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

FourierBodyForce::FourierBodyForce(const fourier_series_data* fourier,
                                   const bool use_circ_model, 
                                   CirculationModel *circ_model,
                                   const INSHierarchyIntegrator* fluid_solver,
                                   const Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
  : d_fourier(fourier), 
    d_use_circ_model(use_circ_model),
    d_circ_model(circ_model),
    d_flux_z(0.0),
    d_fluid_solver(fluid_solver), 
    d_patch_hierarchy(patch_hierarchy)
{

    if (   (( d_use_circ_model) && (!circ_model))
        || ((!d_use_circ_model) && ( circ_model)) ) 
    {
        pout << "Must specify flag for circ model and provide a valid model, or neither.\n"; 
        SAMRAI_MPI::abort();
    }
    
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
    
    // stuff about the physical box and mesh structures
    // take as needed
    const Box<NDIM>& patch_box = patch->getBox();
    
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    
    // Added these
    const double* const x_lower_global = grid_geometry->getXLower();
    const double* const x_upper_global = grid_geometry->getXUpper();
    
    //std::cout << "code thinks the lower domain is " << x_lower_global[0] << ", " << x_lower_global[1] << ", " << x_lower_global[2] << "\n" ;
    //std::cout << "code thinks the upper domain is " << x_upper_global[0] << ", " << x_upper_global[1] << ", " << x_upper_global[2] << "\n" ;
    
    const double z_domain_length    = x_upper_global[2] - x_lower_global[2];

    const double cross_section_area = (x_upper_global[0] - x_lower_global[0]) * 
                                      (x_upper_global[1] - x_lower_global[1]);

    // index without periodicity
    unsigned int k = (unsigned int) floor(data_time / (d_fourier->dt));
    
    // take periodic reduction                         
    unsigned int idx = k % (d_fourier->N_times);

    double force; 

    if (d_use_circ_model){
        // left atrium pressure from circ model 
        const double P_LA = MMHG_TO_CGS * d_circ_model->d_psrc[0];
        
        // left ventricle pressure prescribed from Fourier series 
        const double P_LV = MMHG_TO_CGS * d_fourier->values[idx]; 
        
        // sign is because forward pressure is down 
        force = -(P_LA - P_LV) / z_domain_length; 
    }
    else{ 
        // sign is because forward pressure is down 
        force = -MMHG_TO_CGS * d_fourier->values[idx] / z_domain_length;
    } 

    // Always force in the negative z direction
    const int component = 2;
    //for (int component = 0; component < NDIM; ++component)
    
    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, component)); b; b++){
        const Index<NDIM>& i = b();
        const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);

        (*F_data)(i_s) = force; // FORCE GOES HERE
    }
    
    
    
    // Flow straightener (if desired)
    #ifdef FLOW_STRAIGHTENER
    
        

        // no straightener at time zero
        if (!initial_time){
        
            // some stuff about the geometry
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const x_lower = pgeom->getXLower();
            // const double* const x_upper = pgeom->getXUpper();
    
            // always looking at the z axis here
            const int axis = 2;
    
            // pull velocity variables
            Pointer<SideData<NDIM, double> > U_current_data =
            patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
            Pointer<SideData<NDIM, double> > U_new_data =
            patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
            
            #if !defined(NDEBUG)
                TBOX_ASSERT(U_current_data);
            #endif
            
            // physical height of region of stabilization
            // stabilization is smoothed out from bottom of domain to here
            const double height_physical = 1.0;
            const double max_height_force_applied = height_physical + x_lower_global[2];
            const double center = x_lower_global[axis] + 0.5*height_physical;
    
            // multiplies width of damping by this much
            // full damping support is L*dx
            // required that L*dx + center < height_physical for a continuous mask 
            double L = 1.0; 


            const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
            const double dt = d_fluid_solver->getCurrentTimeStepSize();
            const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
            
            bool full_clamp_on = false; 

            #ifdef FULL_FLOW_CLAMP
                if (data_time < FULL_FLOW_CLAMP_TIME){
                    // F_data->fillAll(0.0);

                    full_clamp_on = true; 

                    // linear decrease in coefficient value 
                    // from max 
                    double k_full_clamp; 
                    if (data_time > 0.0){
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

            if (!full_clamp_on){
                // this may be very, very large
                // consider changing it
                double kappa[NDIM];
                kappa[0] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
                kappa[1] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
                kappa[2] = cycle_num >= 0 ? 0.25 * rho / dt : 0.0; 

                double goal[NDIM]; 
                goal[0] = 0.0; 
                goal[1] = 0.0; 
                
                // z component goal is the average value 
                goal[2] = d_flux_z / cross_section_area; 
                
                // Clamp the velocity in the x,y components
                // Maybe clamp the velocity in the z component too... re
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

                            const double weight    = smooth_kernel((z - center) / (L*dx[axis]*height_physical));
                        
                            (*F_data)(i_s)        += weight*(-kappa[component] * (U - goal[component]));
                            
                            // std::cout << "Placing a force of " << weight*(-kappa * U) << " in component " << component << " at height " << z << "\n";
                            
                        }
                    }
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
