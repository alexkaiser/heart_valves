// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.


// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>
#include <timing.h>
#include <boundary_condition_util.h>
#include <CirculationModel.h>
#include <CirculationModel_with_lv.h>
#include <CirculationModel_aorta.h>
#include <FeedbackForcer.h>
#include <FourierBodyForce.h>

#define NDIM 3 

// #define DEBUG_ELASTICITY 

#define USE_CIRC_MODEL_AORTA

#define MMHG_TO_CGS 1333.22368
#define CGS_TO_MMHG 0.000750061683
#define MPa_TO_CGS 1.0e7


// Elasticity model data.
namespace ModelData
{

    // helper function that is much faster than std::pow(a, -2.0/3.0). std::pow
    // has trouble with some arguments which results in calling slowpow. Since we
    // only call this function on values near 1 (Jacobians) we can be a bit more
    // efficient and compute a*a first since we know that it won't overflow.
    inline
    double fast_pow_n23(const double a)
    {
        return 1.0/std::cbrt(a*a);
    }

    void
    PK1_dev_stress_function(TensorValue<double>& PP,
                            const TensorValue<double>& FF,
                            const libMesh::Point& /*X*/,
                            const libMesh::Point& /*s*/,
                            Elem* const /*elem*/,
                            const vector<const vector<double>*>& /*var_data*/,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* /*ctx*/)
    {

        const TensorValue<double> CC = FF.transpose() * FF;
        const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
        const double J = FF.det();
        const double I1 = CC.tr();

        double mu_s = 1e6; 

        PP.zero();
        PP = mu_s * fast_pow_n23(J) * (FF - (1.0 / 3.0) * I1 * FF_inv_trans);

        return;
    } // PK1_dev_stress_function

    void
    PK1_dil_stress_function(TensorValue<double>& PP,
                            const TensorValue<double>& FF,
                            const libMesh::Point& /*X*/,
                            const libMesh::Point& /*s*/,
                            Elem* const /*elem*/,
                            const vector<const vector<double>*>& /*var_data*/,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double /*time*/,
                            void* /*ctx*/)
    {

        const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
        const double J = FF.det();

        PP.zero();

        double beta_s = 1e6;

        // W(J) = beta_s*(J * log(J) - J + 1)
        PP += beta_s * J * log(J) * FF_inv_trans;

        return;
    } // PK1_dil_stress_function

    // Tether (penalty) force function
    double kappa_s = 1.0e5;
    void tether_force_function(VectorValue<double>& F,
                                const TensorValue<double>& /*FF*/,
                                const libMesh::Point& X,
                                const libMesh::Point& s,
                                Elem* const /*elem*/,
                                const vector<const vector<double>*>& /*var_data*/,
                                const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                double /*time*/,
                                void* /*ctx*/)
    {
        F = kappa_s * (s - X);
        return;
    } // tether_force_function
}
using namespace ModelData;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 MeshBase& mesh,
                 EquationSystems* equation_systems,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/

int main(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // time step used in various places throughout 
        double dt; 
        
        // make some timers
        timestamp_type time1_total, time2_total;    // For total time
        timestamp_type time1, time2;                // For step time
        get_timestamp(&time1_total);                // Get initialization time too

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // check if have a restarted run 
        string restart_read_dirname;
        // int restore_num = 0;
        bool from_restart = false;
        if (argc >= 4)
        {
            // Check whether this appears to be a restarted run.
            FILE* fstream = (SAMRAI_MPI::getRank() == 0 ? fopen(argv[2], "r") : NULL);
            if (SAMRAI_MPI::bcast(fstream != NULL ? 1 : 0, 0) == 1)
            {
                restart_read_dirname = argv[2];
                // restore_num = atoi(argv[3]);
                from_restart = true;
            }
            if (fstream != NULL)
            {
                fclose(fstream);
            }
        }

        int n_restarts_written = 0;
        int max_restart_to_write = 20;
        if (input_db->keyExists("MAX_RESTART_TO_WRITE")){
            max_restart_to_write = input_db->getInteger("MAX_RESTART_TO_WRITE");
        }


        // timestamp_type time1, time2;                // For step time
        double step_time; 

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_filename = app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // std::vector<std::unique_ptr<libMesh::MeshBase> > mesh_vector;
        std::vector<libMesh::MeshBase*> mesh_ptrs;
        std::vector<std::string> part_names;

        Mesh mesh_vessel(init.comm(), NDIM);
        mesh_vessel.read(input_db->getString("MESH_VESSEL"));
        mesh_vessel.prepare_for_use();
        mesh_ptrs.emplace_back(&mesh_vessel);
        part_names.emplace_back(input_db->getString("MESH_VESSEL"));

        Mesh mesh_scaffold(init.comm(), NDIM);
        mesh_scaffold.read(input_db->getString("MESH_SCAFFOLD"));
        mesh_scaffold.prepare_for_use();
        mesh_ptrs.emplace_back(&mesh_scaffold);
        part_names.emplace_back(input_db->getString("MESH_SCAFFOLD"));

        Mesh mesh_valve(init.comm(), NDIM);
        mesh_valve.read(input_db->getString("MESH_VALVE"));
        mesh_valve.prepare_for_use();
        mesh_ptrs.emplace_back(&mesh_valve);
        part_names.emplace_back(input_db->getString("MESH_VALVE"));

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           mesh_ptrs,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.
        IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);

        IBFEMethod::PK1StressFcnData PK1_dev_stress_data_scaffold(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data_scaffold(PK1_dil_stress_function);

        IBFEMethod::PK1StressFcnData PK1_dev_stress_data_valve(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data_valve(PK1_dil_stress_function);

        // added tether 
        IBFEMethod::LagBodyForceFcnData tether_force_data(tether_force_function);
        ib_method_ops->registerLagBodyForceFunction(tether_force_data, 0);

        IBFEMethod::LagBodyForceFcnData tether_force_data_scaffold(tether_force_function);
        ib_method_ops->registerLagBodyForceFunction(tether_force_data_scaffold, 1);

        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
        
        // added part here 
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 0);
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 0);

        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_scaffold, 1);
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_scaffold, 1);

        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data_valve, 2);
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data_valve, 2);

        // added tether 
        ib_method_ops->registerLagBodyForceFunction(tether_force_data, 0);
        ib_method_ops->registerLagBodyForceFunction(tether_force_data_scaffold, 1);

        if (input_db->getBoolWithDefault("ELIMINATE_PRESSURE_JUMPS", false))
        {
            ib_method_ops->registerStressNormalizationPart();
        }
        ib_method_ops->initializeFEEquationSystems();
        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Set up post processor to recover computed stresses.
        ib_method_ops->initializeFEEquationSystems();
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();

        Pointer<IBFEPostProcessor> ib_post_processor =
            new IBFECentroidPostProcessor("IBFEPostProcessor", fe_data_manager);

        ib_post_processor->registerTensorVariable("FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

        pair<IBTK::TensorMeshFcnPtr, void*> PK1_dev_stress_fcn_data(PK1_dev_stress_function, static_cast<void*>(NULL));
        ib_post_processor->registerTensorVariable("sigma_dev",
                                                  MONOMIAL,
                                                  CONSTANT,
                                                  IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                                  vector<SystemData>(),
                                                  &PK1_dev_stress_fcn_data);

        pair<IBTK::TensorMeshFcnPtr, void*> PK1_dil_stress_fcn_data(PK1_dil_stress_function, static_cast<void*>(NULL));
        ib_post_processor->registerTensorVariable("sigma_dil",
                                                  MONOMIAL,
                                                  CONSTANT,
                                                  IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                                                  vector<SystemData>(),
                                                  &PK1_dil_stress_fcn_data);

        Pointer<hier::Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = navier_stokes_integrator->getCurrentContext();
        HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
            /*data_idx*/ -1, "LINEAR_REFINE", /*use_cf_bdry_interpolation*/ false, "CONSERVATIVE_COARSEN", "LINEAR");
        FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                                                QGAUSS,
                                                FIFTH,
                                                /*use_adaptive_quadrature*/ false,
                                                /*point_density*/ 2.0,
                                                /*use_consistent_mass_matrix*/ true,
                                                /*use_nodal_quadrature*/ false);
        ib_post_processor->registerInterpolatedScalarEulerianVariable(
            "p_f", LAGRANGE, FIRST, p_var, p_current_ctx, p_ghostfill, p_interp_spec);

        pout << "to boundary and initial condition creation\n" ; 

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // This is needed to pull variables later
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        #ifndef USE_CIRC_MODEL_AORTA
            // Create Eulerian boundary condition specification objects (when necessary).
            vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));        
            const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
            if (!periodic_domain)
            {   
                pout << "Using b.c. from file, no series for boundary conditions (body force may still have series).\n";
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    ostringstream bc_coefs_name_stream;
                    bc_coefs_name_stream << "u_bc_coefs_" << d;
                    const string bc_coefs_name = bc_coefs_name_stream.str();
                    ostringstream bc_coefs_db_name_stream;
                    bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                    const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                    u_bc_coefs[d] = new muParserRobinBcCoefs(
                        bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
                }
                navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
             
                if (solver_type == "STAGGERED" && input_db->keyExists("BoundaryStabilization"))
                {
                    time_integrator->registerBodyForceFunction(new StaggeredStokesOpenBoundaryStabilizer(
                        "BoundaryStabilization",
                        app_initializer->getComponentDatabase("BoundaryStabilization"),
                        navier_stokes_integrator,
                        grid_geometry));
                }
            }
        #endif


        #ifdef USE_CIRC_MODEL_AORTA

            bool damping_outside = true; 
             
            dt = input_db->getDouble("DT");

            // End systolic / beginning diastolic PA pressure
            double P_aorta_0;
            if (input_db->keyExists("P_aorta_0_MMHG")){
                P_aorta_0 = input_db->getDouble("P_aorta_0_MMHG") * MMHG_TO_CGS;
            }
            else{
                P_aorta_0 = 100.0 * MMHG_TO_CGS;
            }

            std::string fourier_coeffs_name_rv; 
            if (input_db->keyExists("FOURIER_COEFFS_FILENAME_VENTRICLE")){
                fourier_coeffs_name_rv = input_db->getString("FOURIER_COEFFS_FILENAME_VENTRICLE");
            }
            else {
                fourier_coeffs_name_rv = "fourier_coeffs_ventricle.txt"; 
            }
            pout << "attempting to build series with from file: " << fourier_coeffs_name_rv << "\n"; 
            fourier_series_data *fourier_series_ventricle = new fourier_series_data(fourier_coeffs_name_rv.c_str(), dt);

            // boundary vertex files 
            std::string ventricle_vertices_file_name; 
            if (input_db->keyExists("BOUNDARY_FILENAME_VENTRICLE")){
                ventricle_vertices_file_name = input_db->getString("BOUNDARY_FILENAME_VENTRICLE");
            }
            else {
                ventricle_vertices_file_name = "ventricle_bdry.vertex"; 
            }

            std::string aorta_vertices_file_name; 
            if (input_db->keyExists("BOUNDARY_FILENAME_AORTA")){
                aorta_vertices_file_name = input_db->getString("BOUNDARY_FILENAME_AORTA");
            }
            else {
                aorta_vertices_file_name = "aorta_bdry.vertex"; 
            }

            std::string vessel_file_name; 
            if (input_db->keyExists("VESSEL_FILENAME")){
                vessel_file_name = input_db->getString("VESSEL_FILENAME");
            }
            else {
                pout << "VESSEL_FILENAME not found, damping_outside vessel is off\n"; 
                damping_outside = false; 
                vessel_file_name = ""; 
            }

            // scaled cycle length for this patient 
            double t_cycle_length = input_db->getDouble("CYCLE_DURATION");

            // start this far into the Fourier series 
            double t_offset_start_bcs_unscaled; 
            if (input_db->keyExists("T_OFFSET_START_BCS_UNSCALED"))
            {
                t_offset_start_bcs_unscaled = input_db->getDouble("T_OFFSET_START_BCS_UNSCALED"); // starts this 
            }
            else {
                t_offset_start_bcs_unscaled = 0.0; 
            }

            double rcr_on_time; 
            if (input_db->keyExists("RCR_ON_TIME")){
                rcr_on_time = input_db->getDouble("RCR_ON_TIME"); // starts this 
            }
            else {
                rcr_on_time = 0.2; 
            }

            // start in physical time with relation to Fourier series 
            double t_offeset_start = t_offset_start_bcs_unscaled * (t_cycle_length / fourier_series_ventricle->L);

            bool rcr_bcs_on = true;

            bool ventricle_0D_on = input_db->getBoolWithDefault("VENTRICLE_0D_ON", false);

            // start with a linear ramp up in pressure 
            bool P_initial_aorta_equal_to_ventricle = true; 
            // double rcr_on_time = 0.2; 

            CirculationModel_aorta *circ_model_aorta = new CirculationModel_aorta(input_db,
                                                                             fourier_series_ventricle, 
                                                                             ventricle_vertices_file_name,
                                                                             aorta_vertices_file_name,
                                                                             t_cycle_length,
                                                                             t_offset_start_bcs_unscaled, 
                                                                             time_integrator->getIntegratorTime(), 
                                                                             P_aorta_0,
                                                                             rcr_bcs_on, 
                                                                             ventricle_0D_on,
                                                                             P_initial_aorta_equal_to_ventricle, 
                                                                             rcr_on_time); 

            // Create Eulerian boundary condition specification objects.
            vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
            for (int d = 0; d < NDIM; ++d){
                u_bc_coefs[d] = new VelocityBcCoefs_aorta(d, circ_model_aorta);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

            // flow straightener at boundary 
            Pointer<FeedbackForcer> feedback_forcer = new FeedbackForcer(navier_stokes_integrator, patch_hierarchy, NULL, NULL, circ_model_aorta, damping_outside, vessel_file_name, aorta_vertices_file_name);
            time_integrator->registerBodyForceFunction(feedback_forcer);


        #endif // #ifdef USE_CIRC_MODEL_AORTA


        // // Create Eulerian boundary condition specification objects (when necessary).
        // const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        // vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        // if (periodic_shift.min() > 0)
        // {
        //     for (unsigned int d = 0; d < NDIM; ++d)
        //     {
        //         u_bc_coefs[d] = NULL;
        //     }
        // }
        // else
        // {
        //     for (unsigned int d = 0; d < NDIM; ++d)
        //     {
        //         ostringstream bc_coefs_name_stream;
        //         bc_coefs_name_stream << "u_bc_coefs_" << d;
        //         const string bc_coefs_name = bc_coefs_name_stream.str();

        //         ostringstream bc_coefs_db_name_stream;
        //         bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
        //         const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

        //         u_bc_coefs[d] = new muParserRobinBcCoefs(
        //             bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        //     }
        //     navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        // }

        // Create Eulerian body force function specification objects.
        // if (input_db->keyExists("ForcingFunction"))
        // {
        //     Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
        //         "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
        //     time_integrator->registerBodyForceFunction(f_fcn);
        // }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        // std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);
        // set up exodus II IO objects for each part
        std::vector<std::unique_ptr<ExodusII_IO> > exodus_io_vector;
        for(auto mesh_pt : mesh_ptrs)
        {
            exodus_io_vector.emplace_back(uses_exodus ? new ExodusII_IO(*mesh_pt) : nullptr);
        }


        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        if (ib_post_processor) ib_post_processor->initializeFEData();
        pout << "before initializePatchHierarchy\n" ; 
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        pout << "passed initializePatchHierarchy\n" ; 

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                if (ib_post_processor) ib_post_processor->postProcessData(loop_time);
                // exodus_io->write_timestep(
                //     exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                for(int io_idx=0; io_idx<exodus_io_vector.size(); io_idx++)
                {
                    exodus_io_vector[io_idx]->write_timestep(part_names[io_idx], *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
        }

        // // Open streams to save volume of structure.
        // ofstream volume_stream;
        // if (SAMRAI_MPI::getRank() == 0)
        // {
        //     volume_stream.open("volume.curve", ios_base::out | ios_base::trunc);
        // }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        // double dt = 0.0;

        get_timestamp(&time2_total);
        double total_init_time = timestamp_diff_in_seconds(time1_total, time2_total);
        pout << "\n\nTotal initialization time = " << total_init_time << "\n\n"; 
        
        // add some timers         
        get_timestamp(&time1_total); 
        // double step_time; 

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {

            get_timestamp(&time1);          // start step clock 

            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // stop step clock 
            get_timestamp(&time2); 
            step_time = timestamp_diff_in_seconds(time1, time2); 
            
            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "Wallclock time elapsed = " << step_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    if (ib_post_processor) ib_post_processor->postProcessData(loop_time);
                    // exodus_io->write_timestep(
                    //     exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                    for(int io_idx=0; io_idx<exodus_io_vector.size(); io_idx++)
                    {
                        exodus_io_vector[io_idx]->write_timestep(part_names[io_idx], *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                    }
                }
            }

            // Update the circulation model if used 
            #ifdef USE_CIRC_MODEL_AORTA
                {
                    Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                    Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                    Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                    const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                    const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                    Pointer<HierarchyMathOps> hier_math_ops = navier_stokes_integrator->getHierarchyMathOps();
                    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
                    const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
                    circ_model_aorta->advanceTimeDependentData(dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
                }
            #endif

            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            // if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            // {
            //     output_data(patch_hierarchy,
            //                 navier_stokes_integrator,
            //                 mesh,
            //                 equation_systems,
            //                 iteration_num,
            //                 loop_time,
            //                 postproc_data_dump_dirname);
            // }

        } 

        get_timestamp(&time2_total); 
        double total_time = timestamp_diff_in_seconds(time1_total, time2_total); 
        double average_time = total_time / ((double) iteration_num); 
        
        pout << "total run time = " << total_time << " s. \n" ; 
        pout << "average run time = " << average_time << " s. \n" ;
                
        #ifdef USE_CIRC_MODEL_AORTA
            circ_model_aorta->write_plot_code(); 
        #endif 

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            // volume_stream.close();
            std::ofstream done_stream;
            done_stream.open("done.txt", ios_base::out | ios_base::trunc);
            done_stream << "done\n" ;
            done_stream.close(); 
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // run_example

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            MeshBase& mesh,
            EquationSystems* equation_systems,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    file_name = data_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);
    file_name = data_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
} // output_data
