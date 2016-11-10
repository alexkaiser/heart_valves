// Filename: VelocityBcCoefs.C
// Created on 04 May 2007 by Boyce Griffith

// Modified for Fourier series input data 
// 5/2016, Alex Kaiser 

#include "boundary_condition_util.h"
#include "CirculationModel.h"

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
        const Index<NDIM>& i = bc();
        
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
        std::cout << "Failed to open file in Fourier series intialization\n" ; 
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


void fourier_series_data::print_values(){
    // Prints values of the series to command line for debugging 
    
    for (unsigned int j=0; j<N_times; j++) {
        std::cout << values[j] << "\n" ; 
    }    
} 

















