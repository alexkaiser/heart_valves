// Filename: VelocityBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith

// Modified for Fourier series input data 
// 5/2016, Alexander D. Kaiser 

#include "boundary_condition_util.h"
// #include "CirculationModel.h"
#include <CirculationModel_with_lv.h>
#include <CirculationModel_RV_PA.h>
#include <CirculationModel_aorta.h>

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
#include <CartesianPatchGeometry.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>

#include <iostream>
#include <math.h>

#define EXTRA_FWD_PRESSURE

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs(const fourier_series_data *fourier, CirculationModel *circ_model)
    : d_fourier(fourier), d_circ_model(circ_model)
{
    // intentionally blank
    return;
} // VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs

void
VelocityBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                            Pointer<ArrayData<NDIM, double> >& bcoef_data,
                            Pointer<ArrayData<NDIM, double> >& gcoef_data,
                            const Pointer<Variable<NDIM> >& /*variable*/,
                            const Patch<NDIM>& patch,
                            const BoundaryBox<NDIM>& bdry_box,
                            double fill_time) const {
    
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;
    
    // std::cout << "location_index = " << location_index << "\n"; 
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data.isNull());
    #endif
        const Box<NDIM>& bc_coef_box = acoef_data->getBox();
    #if !defined(NDEBUG)
        TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
    #endif

    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const SAMRAI::hier::Index<NDIM>& i = bc();
        
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        if (axis < 2){
            // no slip on x,y axis 
            // std::cout << "no slip on location " << location_index << "\n"; 
            
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
        else if (side == 0){
            
            // Fourier data here
            
            // index without periodicity 
            unsigned int k = (unsigned int) floor(fill_time / (d_fourier->dt));
            
            // take periodic reduction                         
            unsigned int idx = k % (d_fourier->N_times);
            
            a = 0.0; 
            b = 1.0;
            
            // sign for negative in stress tensor
            g = -MMHG_TO_CGS * d_fourier->values[idx];
        
            #ifdef EXTRA_FWD_PRESSURE
                const double extra_fwd_pressure_mmHg = 8.0;
                g += MMHG_TO_CGS * extra_fwd_pressure_mmHg;
            #endif

            //std::cout << "fourier pressure data on location " << location_index << " with value " << g << " or " << fourier->values[idx] << " mmHg\n";

        }
        else if (side == 1){

            a = 0.0;
            b = 1.0;
            // Atrial side pressure from circulation model
            // Stored in CGS units there
            g = -d_circ_model->d_psrc[0];
        }
    }
     
    return;
} // setBcCoefs


IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable
 





/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs_lv_aorta::VelocityBcCoefs_lv_aorta(const int comp_idx,
                                                   CirculationModel_with_lv* circ_model_with_lv)
    : d_comp_idx         (comp_idx), 
      d_circ_model_with_lv(circ_model_with_lv)
{
    // intentionally blank
    return;
} // VelocityBcCoefs_lv_aorta


VelocityBcCoefs_lv_aorta::~VelocityBcCoefs_lv_aorta()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs_lv_aorta


void VelocityBcCoefs_lv_aorta::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                          Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                          Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                          const Pointer<Variable<NDIM> >& /*variable*/,
                                          const Patch<NDIM>& patch,
                                          const BoundaryBox<NDIM>& bdry_box,
                                          double fill_time) const {
    
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;

    //static double last_debug_out = 0.0;  
    
    // std::cout << "location_index = " << location_index << "\n"; 
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data.isNull());
    #endif
        const Box<NDIM>& bc_coef_box = acoef_data->getBox();
    #if !defined(NDEBUG)
        TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
    #endif

    const Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const SAMRAI::hier::Index<NDIM>& i = bc();
        
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        
        if ((axis == 2) && (side == 1) && (d_comp_idx==2)){
            
            // z axis top has all possible interesting boundary conditions 
            // and they all apply to the z component, which has comp_idx 2

            // // always use a time in current cycle 
            // double t_reduced = fill_time - d_cycle_duration * floor(fill_time/d_cycle_duration); 

            // // fourier series has its own period, scale to that 
            // double t_scaled = t_reduced * (d_fourier_aorta->L  / d_cycle_duration); 

            // // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
            // double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

            // // Fourier data here
            // // index without periodicity 
            // unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_aorta->dt));
            
            // // take periodic reduction
            unsigned int idx = d_circ_model_with_lv->d_current_idx_series; 

            /*if (last_debug_out < fill_time){
                pout << "in bc data: t = " << fill_time << ", idx = " << idx << "\n"; 
                d_circ_model_with_lv->print_summary(); 
                last_debug_out = fill_time; 
            }*/ 

            /*
            if (last_debug_out < fill_time){
                pout << "t_scaled_offset = " << t_scaled_offset << "\n"; 
                std::cout << "idx (in Fourier values array) = " << idx << "of " << d_fourier_aorta->N_times <<  "\n";
                std::cout << "MMHG_TO_CGS * d_fourier_aorta->values[idx] = " << MMHG_TO_CGS * d_fourier_aorta->values[idx] << "\n"; 
                std::cout << "MMHG_TO_CGS * d_fourier_atrium->values[idx] = " << MMHG_TO_CGS * d_fourier_atrium->values[idx] << "\n";                         
                last_debug_out = fill_time; 
            } */ 

            double X[NDIM];
            double dist_sq_aorta = 0.0;
            double dist_sq_atrium = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
            }

            const int in_aorta  = d_circ_model_with_lv->point_in_aorta(X[0],X[1]); 
            const int in_atrium = d_circ_model_with_lv->point_in_atrium(X[0],X[1]); 

            if (in_aorta && in_atrium){
                TBOX_ERROR("Position is within both aorta and atrium, should be impossible\n"); 
            }

            if (in_aorta){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -MMHG_TO_CGS * d_circ_model_with_lv->d_fourier_aorta->values[idx];
            }
            else if (in_atrium){
            
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -MMHG_TO_CGS * d_circ_model_with_lv->d_fourier_atrium->values[idx];

            }
            else{
                // no slip outside inlets and outlets 
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }

        }
        else{
            
            // not the z component of the top box 
            // zero tangential slip 
            // zero normal traction 
            if (axis == d_comp_idx){
                // no normal traction for zero pressure 
                a = 0.0;
                b = 1.0;
                g = 0.0;
            }
            else{
                // no tangential slip 
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }

        }
    }
     
    return;
} // setBcCoefs


IntVector<NDIM>
VelocityBcCoefs_lv_aorta::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable



/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs_RV_PA::VelocityBcCoefs_RV_PA(const int comp_idx,
                                                   CirculationModel_RV_PA* circ_model_rv_pa)
    : d_comp_idx        (comp_idx), 
      d_circ_model_rv_pa(circ_model_rv_pa)
{
    // intentionally blank
    return;
} // VelocityBcCoefs_RV_PA


VelocityBcCoefs_RV_PA::~VelocityBcCoefs_RV_PA()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs_RV_PA


void VelocityBcCoefs_RV_PA::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                       Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                       Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                       const Pointer<Variable<NDIM> >& /*variable*/,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double fill_time) const {
    
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;

    //static double last_debug_out = 0.0;  
    
    // std::cout << "location_index = " << location_index << "\n"; 
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data.isNull());
    #endif
        const Box<NDIM>& bc_coef_box = acoef_data->getBox();
    #if !defined(NDEBUG)
        TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
    #endif

    const Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const SAMRAI::hier::Index<NDIM>& i = bc();
        
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        
        if (axis == d_comp_idx){
            // normal component has axis equal to comp_idx 

            // take periodic reduction
            unsigned int idx = d_circ_model_rv_pa->d_current_idx_series; 

            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
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
            else
            {
                TBOX_ERROR("Axis has value that is not 0 1 2\n"); 
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

            if (in_right_ventricle){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -MMHG_TO_CGS * d_circ_model_rv_pa->d_fourier_right_ventricle->values[idx];
                // pout << "Applying pressure of " << d_circ_model_rv_pa->d_fourier_right_ventricle->values[idx] 
                //      << "mmHg to right ventricle at position (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
            }
            else if (in_right_pa){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -d_circ_model_rv_pa->d_right_pa_P; 
                // g = -MMHG_TO_CGS * d_circ_model_rv_pa->d_fourier_right_pa->values[idx];

                // pout << "Applying pressure of " << d_circ_model_rv_pa->d_fourier_right_pa->values[idx] 
                //      << "mmHg to right pa        at position (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
            }
            else if (in_left_pa){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -d_circ_model_rv_pa->d_left_pa_P; 
                // g = -MMHG_TO_CGS * d_circ_model_rv_pa->d_fourier_left_pa->values[idx];
                
                // pout << "Applying pressure of " << d_circ_model_rv_pa->d_fourier_left_pa->values[idx] 
                //      << "mmHg to left pa         at position (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
            }
            else if (side == 0){
                // no slip on lower boundary (all inlets and outlets are lower)
                a = 1.0;
                b = 0.0;
                g = 0.0;                
            }
            else{
                // no normal traction for zero pressure 
                a = 0.0;
                b = 1.0;
                g = 0.0;
            }
        
        }
        else {
            // no tangential slip everywhere on domain 
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }

    }
     
    return;
} // setBcCoefs


IntVector<NDIM>
VelocityBcCoefs_RV_PA::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable





/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs_aorta::VelocityBcCoefs_aorta(const int comp_idx,
                                             CirculationModel_aorta* circ_model_aorta)
    : d_comp_idx        (comp_idx), 
      d_circ_model_aorta(circ_model_aorta)
{
    // intentionally blank
    return;
} // VelocityBcCoefs_aorta


VelocityBcCoefs_aorta::~VelocityBcCoefs_aorta()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs_aorta


void VelocityBcCoefs_aorta::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                       Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                       Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                       const Pointer<Variable<NDIM> >& /*variable*/,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double fill_time) const {
    
    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;

    //static double last_debug_out = 0.0;  
    
    // std::cout << "location_index = " << location_index << "\n"; 
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data.isNull());
    #endif
        const Box<NDIM>& bc_coef_box = acoef_data->getBox();
    #if !defined(NDEBUG)
        TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
    #endif

    const Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const SAMRAI::hier::Index<NDIM>& i = bc();
        
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        
        if (axis == d_comp_idx){
            // normal component has axis equal to comp_idx 

            // take periodic reduction
            unsigned int idx = d_circ_model_aorta->d_current_idx_series; 

            double X[NDIM];
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
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
            else
            {
                TBOX_ERROR("Axis has value that is not 0 1 2\n"); 
            }

            const int in_ventricle  = d_circ_model_aorta->point_in_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
            const int in_aorta      = d_circ_model_aorta->point_in_aorta    (X_in_plane_1, X_in_plane_2, axis, side);

            if (in_ventricle && in_aorta){
                TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
            }

            if (in_ventricle){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                g = -d_circ_model_aorta->d_ventricle_P; 
                // pout << "Applying pressure of " << d_circ_model_rv_pa->d_fourier_right_ventricle->values[idx] 
                //      << "mmHg to right ventricle at position (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
            }
            else if (in_aorta){
                // sign for negative in stress tensor
                a = 0.0; 
                b = 1.0; 
                if (d_circ_model_aorta->d_rcr_bcs_on){
                    g = -d_circ_model_aorta->d_aorta_P; 
                }
                else{
                    TBOX_ERROR("not implemented\n");
                }
                // pout << "Applying pressure of " << d_circ_model_rv_pa->d_fourier_right_pa->values[idx] 
                //      << "mmHg to right pa        at position (" << X[0] << ", " << X[1] << ", " << X[2] << ")\n"; 
            }
            else if ((side == d_circ_model_aorta->d_aorta_side) && (axis == d_circ_model_aorta->d_aorta_axis)){
                // no slip on relevant boundary
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }
            else if ((side == d_circ_model_aorta->d_ventricle_side) && (axis == d_circ_model_aorta->d_ventricle_axis)){
                // no slip on relevant boundary
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }
            else{
                // no normal traction for zero pressure 
                a = 0.0;
                b = 1.0;
                g = 0.0;
            }
        
        }
        else {
            // no tangential slip everywhere on domain 
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }

    }
     
    return;
} // setBcCoefs


IntVector<NDIM>
VelocityBcCoefs_aorta::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable




fourier_series_data::fourier_series_data(string file_name, double dt_input){
    /*
    Constructor 
    
    Reads Fourier series data from file. 
    
    Input: 
        file_name   Text file with Fourier coefficients 
        dt          Simulation time step 

    File format: 
        N           Number of Fourier coefficients 
        L           Period of Fourier series 
        a_0         Cosine zero coeff
        a_1 b_1     Cosine and sine coeffs, must be N
        ...
    */ 
    
    ifstream f; 
    f.open(file_name.c_str(), ios::in);
    
    if (!f.is_open()){
        pout << "file_name.c_str() = " << file_name.c_str() << "\n"; 
        TBOX_ERROR("Failed to open file in Fourier series intialization\n"); 
    } 
    
    // set the local dt 
    dt = dt_input; 
    
    // number of coefficients 
    unsigned int N; 
    f >> N; 
    
    // Length of interval 
    f >> L; 
    
    double *a_n = new double[N]; 
    double *b_n = new double[N]; 
    
    double pi = 4.0*atan2(1,1); 
        
    f >> a_n[0]; 
    b_n[0] = 0.0; 
    for(unsigned int j=1; (j<N) && (!f.eof()); j++){
        f >> a_n[j]; 
        f >> b_n[j]; 
    }
    
    N_times = (unsigned int) floor(L/dt);
     
    values = new double[N_times]; 
    
    unsigned int t_idx; 
    double t;
    
    for (t_idx=0, t=0.0; t_idx < N_times; t_idx++, t+=dt){
        for (unsigned int j=0; j<N; j++){
            values[t_idx] += a_n[j] * cos((2*pi/L) * j * t); 
            values[t_idx] += b_n[j] * sin((2*pi/L) * j * t); 
        }     
    }
    
    f.close();     
    delete[] a_n; 
    delete[] b_n;      
} 

fourier_series_data::~fourier_series_data(){
    // Basic destructor, deletes values array 
    if (values)
        delete[] values; 
} 


void fourier_series_data::print_values() const{
    // Prints values of the series to command line for debugging 
    
    pout << "N_times = " << N_times << "\n"; 
    pout << "L = " << L << "\n"; 
    pout << "dt = " << dt << "\n"; 
    for (unsigned int j=0; j<N_times; j++) {
        pout << values[j] << "\n" ; 
    }    
} 



ventricle_0D_model::ventricle_0D_model(Pointer<Database> input_db, 
                                       double cycle_duration, 
                                       double initialization_time)
    :
    d_object_name("ventricle_0D_model_name"),  // constant name here  
    d_registered_for_restart(true),      // always true
    d_cycle_duration(cycle_duration),
    d_initialization_time(initialization_time)
    
    {

    double dt = input_db->getDouble("DT");

    std::string fourier_coeffs_name_Q_in = input_db->getStringWithDefault("FOURIER_COEFFS_Q_IN" , "fourier_coeffs_Q_mi.txt"); 
    std::string fourier_coeffs_name_act  = input_db->getStringWithDefault("FOURIER_COEFFS_ACT" , "fourier_coeffs_lv_activation.txt"); 

    d_fourier_q_in_ventricle = new fourier_series_data(fourier_coeffs_name_Q_in.c_str(), dt);
    d_act_ventricle = new fourier_series_data(fourier_coeffs_name_act, dt); 

    d_V_rest_diastole = input_db->getDoubleWithDefault("V_REST_DIASTOLE", 26.1);
    d_V_rest_systole  = input_db->getDoubleWithDefault("V_REST_SYSTOLE", 18.0);
    d_E_min           = input_db->getDoubleWithDefault("E_MIN", 1.818181818181818e+02);
    d_E_max           = input_db->getDoubleWithDefault("E_MAX", 1.057082452431290e+04); 
    d_R_LVOT          = input_db->getDoubleWithDefault("R_LVOT", 0.0); 
    d_inductance      = input_db->getDoubleWithDefault("INDUCTANCE", 0.0); 

    double V_initial = input_db->getDoubleWithDefault("V_INITIAL", d_V_rest_diastole);

    const bool print_debug = true; 
    if (print_debug){
        pout << "ventricle_0D_model initializing.\n"; 
        pout << "V_initial = " << V_initial << "\n";
        pout << "d_cycle_duration = " << d_cycle_duration << "\n"; 
        pout << "d_initialization_time = " << d_initialization_time << "\n"; 
        pout << "d_V_rest_diastole = " << d_V_rest_diastole  << "\n"; 
        pout << "d_V_rest_systole = " << d_V_rest_systole << "\n"; 
        pout << "d_E_min = " << d_E_min << "\n"; 
        pout << "d_E_max = " << d_E_max << "\n"; 
        pout << "d_R_LVOT = " << d_R_LVOT << "\n";
        pout << "d_inductance = " << d_inductance << "\n"; 
        //pout << "d_fourier_q_in_ventricle = " << "\n"; 
        // d_fourier_q_in_ventricle->print_values(); 
        //pout << "d_act_ventricle = " << "\n"; 
        //d_act_ventricle->print_values(); 
    }


    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }    
    else {
        // default initialization not from restart 
        d_Q_out = 0.0; 
        d_Q_out_prev = 0.0; 
        d_time = 0.0; 
        d_current_idx_series = 0; 

        double d_act_temp = 0.0; 

        // use equations with zero activation  
        d_V_ventricle = V_initial; 

        // rest volume  
        d_V_rest_ventricle = (1.0 - d_act_temp) * (d_V_rest_diastole - d_V_rest_systole) + d_V_rest_systole; 

        // elastance 
        d_Elas = (d_E_max - d_E_min) * d_act_temp + d_E_min; 

        // pressure 
        d_P_ventricle_in = d_Elas * (d_V_ventricle - d_V_rest_ventricle); 
        d_P_lvot_upstream = d_P_ventricle_in;  
        d_P_ventricle = d_P_lvot_upstream;  
    }

}

ventricle_0D_model::~ventricle_0D_model(){

    if(d_fourier_q_in_ventricle) {
        delete d_fourier_q_in_ventricle; 
    }
    if(d_act_ventricle){
        delete d_act_ventricle; 
    }

}


void ventricle_0D_model::advanceTimeDependentData(double dt, double time, double Q_out){

    /* 
    % eqn 3 for ventricular volume 
    V_lv(n+1) = V_lv(n) + dt * (Q_mi(n) - Q_ao(n)); 

    % elastance updates 
    Vrest(n+1) = (1.0 - act(n+1)) * (Vrd - Vrs) + Vrs;
    Elas(n+1) = (Emax - Emin) * act(n+1) + Emin;

    % pressure update 
    P_lv_in(n+1) = Elas(n+1) * (V_lv(n+1) - Vrest(n+1)); 
    */ 

    // save last time step flow 
    d_Q_out_prev = d_Q_out; 

    d_Q_out = Q_out; 

    d_time = time; 

    // compute which index in the Fourier series we need here 
    // always use a time in current cycle 
    double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration); 

    // fourier series has its own period, scale to that 
    double t_scaled = t_reduced * (d_fourier_q_in_ventricle->L  / d_cycle_duration); 

    // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
    // no scaling allowed for now 
    double t_scaled_offset = t_scaled; // + d_t_offset_bcs_unscaled; 

    // Fourier data here
    // index without periodicity 
    unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_q_in_ventricle->dt));
    
    // // take periodic reduction
    d_current_idx_series = k % (d_fourier_q_in_ventricle->N_times);

    d_Q_in = d_fourier_q_in_ventricle->values[d_current_idx_series]; 

    // volume update 
    d_V_ventricle += dt * (d_Q_in - d_Q_out); 

    d_act_temp = 0.0; 
    if (d_time > d_initialization_time){
        d_act_temp = d_act_ventricle->values[d_current_idx_series];
    }

    // rest volume  
    d_V_rest_ventricle = (1.0 - d_act_temp) * (d_V_rest_diastole - d_V_rest_systole) + d_V_rest_systole;

    // elastance 
    d_Elas = (d_E_max - d_E_min) * d_act_temp + d_E_min;

    // pressure 
    d_P_ventricle_in = d_Elas * (d_V_ventricle - d_V_rest_ventricle);

    // inductor 
    if (d_time > d_initialization_time){
        d_P_lvot_upstream = d_P_ventricle_in - d_inductance * (d_Q_out - d_Q_out_prev)/dt;
    }
    else{
        double inductance_tmp = (d_time / d_initialization_time) * d_inductance; 
        d_P_lvot_upstream = d_P_ventricle_in - inductance_tmp * (d_Q_out - d_Q_out_prev)/dt;
    }

    d_P_ventricle = d_P_lvot_upstream - d_R_LVOT * d_Q_out; 

}


void ventricle_0D_model::putToDatabase(Pointer<Database> db)
{

    db->putDouble("d_Q_out", d_Q_out);
    db->putDouble("d_Q_out_prev", d_Q_out_prev);
    db->putDouble("d_V_ventricle", d_V_ventricle);
    db->putDouble("d_V_rest_ventricle", d_V_rest_ventricle);
    db->putDouble("d_P_ventricle_in", d_P_ventricle_in);
    db->putDouble("d_P_lvot_upstream", d_P_lvot_upstream);
    db->putDouble("d_P_ventricle", d_P_ventricle);

    db->putDouble("d_Elas", d_Elas);

    db->putDouble("d_Q_in", d_Q_in);
    db->putDouble("d_act_temp", d_act_temp);

    db->putDouble("d_initialization_time", d_initialization_time);
    db->putDouble("d_time", d_time);
    db->putDouble("d_cycle_duration", d_cycle_duration);
    db->putInteger("d_current_idx_series", d_current_idx_series);

    return; 
} // putToDatabase


void ventricle_0D_model::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_Q_out               = db->getDouble("d_Q_out"); 
    d_Q_out_prev          = db->getDouble("d_Q_out_prev"); 
    d_V_ventricle         = db->getDouble("d_V_ventricle"); 
    d_V_rest_ventricle    = db->getDouble("d_V_rest_ventricle"); 
    d_P_ventricle_in      = db->getDouble("d_P_ventricle_in"); 
    d_P_lvot_upstream     = db->getDouble("d_P_lvot_upstream");
    d_P_ventricle         = db->getDouble("d_P_ventricle"); 
    d_Elas                = db->getDouble("d_Elas"); 

    d_Q_in                = db->getDouble("d_Q_in");
    d_act_temp            = db->getDouble("d_act_temp");

    d_initialization_time = db->getDouble("d_initialization_time"); 
    d_time                = db->getDouble("d_time");
    d_cycle_duration      = db->getDouble("d_cycle_duration"); 
    d_current_idx_series  = db->getInteger("d_current_idx_series"); 

    return;
} // getFromRestart





