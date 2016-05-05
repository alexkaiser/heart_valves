// Filename: VelocityBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith

// Modified for Fourier series input data 
// 5/2016, Alex Kaiser 

#include "boundary_condition_util.h"

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

#include <iostream>
#include <math.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

VelocityBcCoefs::VelocityBcCoefs(const fourier_series_data *fourier)
    // : can put some default stuff here  
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
                            double /*fill_time*/) const
{

    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;
    
    #if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data.isNull());
    #endif
        const Box<NDIM>& bc_coef_box = acoef_data->getBox();
    #if !defined(NDEBUG)
        TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
    #endif
    
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    //const double* const dx = pgeom->getDx();
    //const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const Index<NDIM>& i = bc();
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        if (axis != 1)
        {
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
     /*   else if (d_comp_idx != axis)
        {
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
        else if (d_comp_idx == axis) // FIX BOUNDARY 
        {

            a = (r <= rsrc) ? 0.0 : 1.0;
            b = (r <= rsrc) ? 1.0 : 0.0;
            g = (r <= rsrc) ? -psrc : 0.0;
        }
        */ 
    }

    
    /*const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const int side = location_index % 2;
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const Index<NDIM>& i = bc();
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);
        if (axis != 1)
        {
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
        else if (d_comp_idx != axis)
        {
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }
        else if (d_comp_idx == axis)
        {
            const double psrc = d_circ_model->d_psrc[side];
            const double rsrc = d_circ_model->d_rsrc[side];
            const Point& posn = d_circ_model->d_posn[side];
            double X[NDIM];
            double r_sq = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                if (d != axis) r_sq += pow(X[d] - posn[d], 2.0);
            }
            const double r = sqrt(r_sq);
            a = (r <= rsrc) ? 0.0 : 1.0;
            b = (r <= rsrc) ? 1.0 : 0.0;
            g = (r <= rsrc) ? -psrc : 0.0;
        }
    }
    */ 
    return;
} // setBcCoefs

IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable




fourier_series_data::fourier_series_data(string file_name, double dt){
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

    std::cout << "in constructor \n" ;     
    
    ifstream f; 
    f.open(file_name.c_str(), ios::in);
    
    if (!f.is_open()){
        std::cout << "Failed to open file in Fourier series intialization\n" ; 
    } 
    
    std::cout << "file is open and if statement passed \n" ;     
    
    // number of coefficients 
    int N; 
    f >> N; 
    
    // Length of interval 
    f >> L; 
    
    double *a_n = new double[N]; 
    double *b_n = new double[N]; 
    
    double pi = 4.0*atan2(1,1); 
        
    std::cout << "to read loop \n" ; 
        
    f >> a_n[0]; 
    b_n[0] = 0.0; 
    for(int j=1; (j<N) && (!f.eof()); j++){
        f >> a_n[j]; 
        f >> b_n[j]; 
    }
    
    N_times = (int) floor(L/dt);
     
    values = new double[N]; 
    
    int t_idx; 
    double t;

    std::cout << "to evel loop \n" ; 
    
    for (t_idx=0, t=0.0; t_idx < N_times; t_idx++, t+=dt){
        for (int j=0; j<N; j++){
            values[t_idx] += a_n[j] * cos((2*pi/L) * j * t); 
            values[t_idx] += b_n[j] * sin((2*pi/L) * j * t); 
        }     
    }

    std::cout << "to cleanup \n" ; 
    
    f.close();     
    delete[] a_n; 
    delete[] b_n; 
     
} 

fourier_series_data::~fourier_series_data(){
    // Basic destructor, deletes values array 
    if (values)
        delete[] values; 
} 


void fourier_series_data::print_values(){
    // Prints values of the series to command line for debugging 
    
    for (int j=0; j<N_times; j++) {
        std::cout << values[j] << "\n" ; 
    }    
} 

















