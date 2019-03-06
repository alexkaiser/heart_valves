// Filename: CirculationModel_with_lv.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_with_lv.h"

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
#include <PatchLevel.h>
#include <SideData.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

#include <Eigen/Dense>
using namespace Eigen;

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel_with_lv::CirculationModel_with_lv(const fourier_series_data *fourier_aorta, 
                                                   const fourier_series_data *fourier_atrium, 
                                                   const double  radius_aorta,
                                                   const double  radius_atrium,
                                                   const double *center_aorta,
                                                   const double *center_atrium, 
                                                   const double cycle_duration,
                                                   const double t_offset_bcs_unscaled)
    : d_fourier_aorta (fourier_aorta), 
      d_fourier_atrium(fourier_atrium), 
      d_radius_aorta  (radius_aorta),
      d_radius_atrium (radius_atrium),
      d_center_aorta  (center_aorta),
      d_center_atrium (center_atrium),       
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_Q_aorta       (0.0), 
      d_Q_left_atrium (0.0)
{
    // intentionally blank
    return;
} // CirculationModel

CirculationModel_with_lv::~CirculationModel_with_lv()
{
    return;
} // ~CirculationModel_with_lv


void CirculationModel_with_lv::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    
    double Q_aorta_local = 0.0; 
    double Q_left_atrium_local = 0.0; 

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }
                double X[NDIM];

                // always looking for z flux here 
                // side is always 1, top of box 
                static const int axis = 2;
                int side = 1; 
                const bool is_lower = (side == 0);
                if (pgeom->getTouchesRegularBoundary(axis, side))
                {

//                    radius and position of sources are handled by input data 
//                    const double rsrc = d_rsrc[side];
//                    const Point& posn = d_posn[side];
                    
                    Vector n;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                    }
                    Box<NDIM> side_box = patch_box;
                    if (is_lower)
                    {
                        side_box.lower(axis) = patch_box.lower(axis);
                        side_box.upper(axis) = patch_box.lower(axis);
                    }
                    else
                    {
                        side_box.lower(axis) = patch_box.upper(axis) + 1;
                        side_box.upper(axis) = patch_box.upper(axis) + 1;
                    }
                    for (Box<NDIM>::Iterator b(side_box); b; b++)
                    {
                        const Index<NDIM>& i = b();

                        double X[NDIM];
                        double dist_sq_aorta = 0.0;
                        double dist_sq_atrium = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                            if (d != axis){
                                dist_sq_aorta  += pow(X[d] - d_center_aorta[d],  2.0);
                                dist_sq_atrium += pow(X[d] - d_center_atrium[d], 2.0);
                            }
                        }
                        const double dist_aorta  = sqrt(dist_sq_aorta);
                        const double dist_atrium = sqrt(dist_sq_atrium);

                        if (dist_aorta < d_radius_aorta)
                        {
                            const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                            if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                            {
                                double dA = n[axis] * dV / dx[axis];
                                Q_aorta_local += (*U_data)(i_s)*dA;
                            }
                        }

                        if (dist_atrium < d_radius_atrium)
                        {
                            const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                            if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                            {
                                double dA = n[axis] * dV / dx[axis];
                                Q_left_atrium_local += (*U_data)(i_s)*dA;
                            }
                        }
                    }
                }
                
            }
        }
    }

    d_Q_aorta       = SAMRAI_MPI::sumReduction(Q_aorta_local);
    d_Q_left_atrium = SAMRAI_MPI::sumReduction(Q_left_atrium_local);

} // advanceTimeDependentData


void
CirculationModel_with_lv::putToDatabase(Pointer<Database> db)
{
    return; 
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CirculationModel_with_lv::writeDataFile() const
{
    return;
} // writeDataFile

void
CirculationModel_with_lv::getFromRestart()
{
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////