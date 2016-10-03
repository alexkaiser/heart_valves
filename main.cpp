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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>



// added this header
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/IBInstrumentPanel.h>
#include <ibamr/IBStandardSourceGen.h>
#include <vector>
#include <queue>
#include <cmath>
#include <timing.h>
#include <boundary_condition_util.h>
#include <CirculationModel.h>
#include <FeedbackForcer.h>
#include <FourierBodyForce.h>


// #define IMPLICIT_SOLVER
#ifdef IMPLICIT_SOLVER
     #include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#endif


#if defined(IBAMR_HAVE_SILO)
#include <silo.h>
#endif

//void init_source_variables(vector<double>& Q_src, const vector<double>& P_src);
//void set_source_variables(vector<double>& Q_src, const vector<double>& P_src, double current_time, double dt, bool sink_on);
//void update_rest_lengths(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* const l_data_manager, const double alpha);

void update_target_point_positions(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* const l_data_manager, const double current_time, const double dt, fourier_series_data *fourier_body_force);

inline double spring_function_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);

#define DEBUG_OUTPUT 0 
#define ENABLE_INSTRUMENTS
#define FOURIER_SERIES_BC
#define FOURIER_SERIES_BODY_FORCE
#define DYNAMIC_BOUNDARY_STAB




#define MMHG_TO_CGS 1333.22368
#define CGS_TO_MMHG 0.000750061683
#define MPa_TO_CGS 1.0e7

#define SOURCE_ON_TIME       0.1
#define CONST_SRC_TIME       0.3
#define CONST_SRC_STRENGTH  93.0

#define MAX_STEP_FOR_CHANGE 1000

// #define MOVING_PAPILLARY


namespace{
    inline double smooth_kernel(const double r){
        return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
    } // smooth_kernel
}



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
int main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
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

        #ifdef IMPLICIT_SOLVER
            // Read default Petsc options
            if (input_db->keyExists("petsc_options_file"))
            {
                std::string PetscOptionsFile = input_db->getString("petsc_options_file");
                #if (!PETSC_VERSION_RELEASE)
                    PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, PetscOptionsFile.c_str(), PETSC_TRUE);
                #else
                    PetscOptionsInsertFile(PETSC_COMM_WORLD, PetscOptionsFile.c_str(), PETSC_TRUE);
                #endif
            }
        #endif

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
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
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        
        
        #ifdef IMPLICIT_SOLVER
            Pointer<IBHierarchyIntegrator> time_integrator =
                new IBImplicitStaggeredHierarchyIntegrator("IBHierarchyIntegrator",
                                                            app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                            ib_method_ops,
                                                            navier_stokes_integrator);
        #else
            Pointer<IBHierarchyIntegrator> time_integrator =
                new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                                    app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                    ib_method_ops,
                                                    navier_stokes_integrator);
        #endif
        
        
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

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        
        // adding custom function for collagen springs
        ib_force_fcn->registerSpringForceFunction(1, &spring_function_collagen, &deriv_spring_collagen);
        
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);
        
        
        #ifdef ENABLE_INSTRUMENTS
         
            std::vector<double> flux_valve_ring;
            std::ofstream flux_output_stream;
            
            // no idea what these do 
            int u_data_idx = 0; 
            int p_data_idx = 0; 
            
            Pointer<IBInstrumentPanel> instruments = new IBInstrumentPanel("meter_0", input_db);
        #endif
        
        
        
        LDataManager* l_data_manager = ib_method_ops->getLDataManager();
        pout << "passed LDataManager creation\n" ; 

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

        // Create Eulerian boundary condition specification objects (when necessary).
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));
        
        #ifdef FOURIER_SERIES_BC
            pout << "to Fourier series creation\n"; 
        
            dt = input_db->getDouble("DT"); 
            
            pout << "to constructor\n"; 
            fourier_series_data *fourier_ventricle = new fourier_series_data("fourier_coeffs_ventricle.txt", dt);
            pout << "series data successfully built\n";
        
            // fourier->print_values();
        
            // Need to declare these out here for scope reasons
            // index without periodicity
            double start = input_db->getDouble("START_TIME");
            unsigned int k = (unsigned int) floor(start / (fourier_ventricle->dt));
        
            // take periodic reduction
            unsigned int idx = k % (fourier_ventricle->N_times);
        
            bool restart_circ_model = true;
    
            // End systolic / beginning diastolic PA pressure
            double P_PA_0 = 18;
    
            // Simulation starts at moment when LA and LV pressures are equal
            // Note that circ model has units of mmHg
            // As does Fourier series 
            double P_LA_0 = fourier_ventricle->values[idx];
    
            // Circulation model
            CirculationModel *circ_model   = new CirculationModel("circ_model", P_PA_0, P_LA_0, start, restart_circ_model);
        #endif
        
        // This is needed to pull variables later
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
        if (!periodic_domain)
        {
        
            #ifdef FOURIER_SERIES_BC
                
                if (NDIM != 3){
                    pout << "Current implementation requires 3d\n"; 
                    SAMRAI_MPI::abort();  
                }
                    
                for (unsigned int d = 0; d < 2; ++d)
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
            
                // manually update third component, 
                // which is the only one not easily set in the input file
                VelocityBcCoefs *z_bdry_coeffs = new VelocityBcCoefs(fourier_ventricle, circ_model);
                u_bc_coefs[2] = z_bdry_coeffs;
            
            
                #ifdef DYNAMIC_BOUNDARY_STAB
            
                    const bool z_periodic = (grid_geometry->getPeriodicShift())[2];
                    if (z_periodic){
                        pout << "Code thinks that z is periodic, outflow boundary stabilization off\n" ;
                    }
            
                    // only for staggered grids
                    if ((solver_type == "STAGGERED") && (!z_periodic)){
                        // always the Z component
                        Pointer<FeedbackForcer> bdry_dynamic_stab = new FeedbackForcer(z_bdry_coeffs, navier_stokes_integrator, patch_hierarchy);
                        time_integrator->registerBodyForceFunction(bdry_dynamic_stab);
                    }
                #endif
            
            #else 
                pout << "Using b.c. from file, no Fourier series created.\n";
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
            #endif 


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

        // generic body force
        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            if (input_db->keyExists("BoundaryStabilization"))
            {
                TBOX_ERROR("Cannot currently use boundary stabilization with additional body forcing");
            }
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }


        #ifdef FOURIER_SERIES_BODY_FORCE
            const bool z_periodic = (grid_geometry->getPeriodicShift())[2];
            if (z_periodic){
                pout << "to Fourier series creation with body force\n";
            
                dt = input_db->getDouble("DT");
            
                pout << "to constructor\n";
                fourier_series_data *fourier_body_force = new fourier_series_data("fourier_coeffs.txt", dt);
                pout << "series data successfully built\n";
        
                Pointer<FourierBodyForce> body_force = new FourierBodyForce(fourier_body_force, navier_stokes_integrator, patch_hierarchy);
                time_integrator->registerBodyForceFunction(body_force);
            }
        #endif


        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        pout << "before initializePatchHierarchy\n" ; 
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        pout << "passed initializePatchHierarchy\n" ; 

        // Finest level does not change throughout 
        const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
        

        #ifdef ENABLE_INSTRUMENTS
            // do this after initialize patch hierarchy
            instruments->initializeHierarchyIndependentData(patch_hierarchy, l_data_manager);
            
            // Stream to write-out flux data 
            if (!from_restart){
                if (SAMRAI_MPI::getRank() == 0){
                    flux_output_stream.open("flux_plot_IBAMR.m", ios_base::out | ios_base::trunc);
                    flux_output_stream << "data = [" ; 
                }
            }
            // if we are restarting, we want to append to the file 
            // we are assuming that the pressure plot write was not interuppted here in the previous run 
            else{
                if (SAMRAI_MPI::getRank() == 0){
                    flux_output_stream.open("flux_plot_IBAMR.m", ios_base::out | ios_base::app);
                }
            }
            
        #endif


        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Setup Silo writer.
        if (silo_data_writer)
        {
            //const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
            Pointer<LData> F_data = l_data_manager->getLData("F", finest_hier_level);
            silo_data_writer->registerVariableData("F", F_data, finest_hier_level);
        }

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        
        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
            
            #ifdef ENABLE_INSTRUMENTS
                
                if (instruments->isInstrumented())
                    pout << "is instrumented passed\n" ;
                else
                    pout << "is instrumented returned FALSE\n" ;  
                    
                instruments->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                instruments->readInstrumentData(u_data_idx, p_data_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                
                flux_valve_ring = instruments->getFlowValues(); 
                // pout << "flux at t = " << loop_time << ", Q = " << flux_valve_ring[0] << "\n"; 
                
                if (SAMRAI_MPI::getRank() == 0){
                    flux_output_stream << loop_time << ",\t" << flux_valve_ring[0] << ";\n"; 
                    flux_output_stream.flush(); 
                }
                
            #endif
            
        }
 

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        dt = 0.0;
        double dt_prev = 0.0;
        bool prev_step_initialized = false;
        
        get_timestamp(&time2_total);
        double total_init_time = timestamp_diff_in_seconds(time1_total, time2_total);
        pout << "\n\nTotal initialization time = " << total_init_time << "\n\n"; 
        
        // add some timers         
        get_timestamp(&time1_total); 
        double step_time; 

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            
            get_timestamp(&time1);          // start step clock 
            
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();


            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            
            // save the last time step 
            dt_prev = dt; 
            
            // get current step 
            dt = time_integrator->getMaximumTimeStepSize();
        
            if((!prev_step_initialized) && from_restart){
                    dt_prev = dt;
                    prev_step_initialized = true;
            }
            
            // update target locations if they are moving
            #ifdef MOVING_PAPILLARY
                #ifdef FOURIER_SERIES_BODY_FORCE
                    update_target_point_positions(patch_hierarchy, l_data_manager, loop_time, dt, fourier_body_force);
                #else
                    pout << "other papillary movement not implemented, must use periodic with fourier series for now\n";
                    SAMRAI_MPI::abort();
                #endif
            #endif
            
            // In non periodic case, update the circulation model
            const bool z_periodic = (grid_geometry->getPeriodicShift())[2]; 
            if (!z_periodic){
                Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                Pointer<HierarchyMathOps> hier_math_ops = navier_stokes_integrator->getHierarchyMathOps();
                const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
                const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
                circ_model->advanceTimeDependentData(
                    dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
            } 
            
            // step the whole thing
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
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
                
            }
            
            // Write this every step
            #ifdef ENABLE_INSTRUMENTS
                
                Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                
                instruments->initializeHierarchyDependentData(patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                instruments->readInstrumentData(U_current_idx, P_current_idx, patch_hierarchy, l_data_manager, iteration_num, loop_time); 
                flux_valve_ring = instruments->getFlowValues(); 
                // pout << "flux at t = " << loop_time << ", Q = " << flux_valve_ring[0] << "\n";
                
                if (SAMRAI_MPI::getRank() == 0){
                    flux_output_stream << loop_time << ",\t" << flux_valve_ring[0] << ";\n"; 
                    flux_output_stream.flush(); 
                }                
            
            #endif
            
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            
            
            if (iteration_num > MAX_STEP_FOR_CHANGE){
                // if there is a change 
                if (dt != dt_prev){
                    // ignore the first and last steps  
                    if (!last_step){
                        pout << "Timestep change encountered (manual, after max change step).\n" ;
                        pout << "dt = " << dt << ",\tdt_prev = " << dt_prev << "\n";
                        SAMRAI_MPI::abort();
                    }
                }
            }
                        

            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            
        }

        get_timestamp(&time2_total); 
        double total_time = timestamp_diff_in_seconds(time1_total, time2_total); 
        double average_time = total_time / ((double) iteration_num); 
        
        pout << "total run time = " << total_time << " s. \n" ; 
        pout << "average run time = " << average_time << " s. \n" ;
        
        #ifdef ENABLE_INSTRUMENTS
            if (SAMRAI_MPI::getRank() == 0){
                flux_output_stream << "]; \n\n"; 
                
                flux_output_stream << "fig = figure;\n subplot(2,1,1);\n plot(data(:,1), -data(:,2));\n xlabel('t');\n ylabel('flux (cm^3 / s)  ');\n";
                flux_output_stream << "hold on;\n plot(data(:,1), 0*data(:,2)); \n"; 
                flux_output_stream << "dt = " << dt << "; \n"; 
                flux_output_stream << "net_flux = dt*cumsum(-data(:,2));\n subplot(2,1,2);\n plot(data(:,1), net_flux);\n xlabel('t');\n ylabel('net flux (cm^3)');\n"; 
                flux_output_stream << "printfig(fig,'flux.eps');\n"; 
                
                flux_output_stream.close();
            }
        #endif

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main



void update_target_point_positions(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* const l_data_manager, const double current_time, const double dt, fourier_series_data *fourier_body_force){
    // Requires to have pressure difference as a Fourier series
    // Atrial pressure is positive if higher
    // so positive pressure drives forward flow

    // quick return at beginning so that things do not move in discontinuous manner
    if (current_time < 0.1)
        return; 


    const static double LEFT_PAPILLARY[3]  = {-0.972055648767080, -1.611924550017006, -2.990100960298683};
    const static double RIGHT_PAPILLARY[3] = {-1.542417595752084,  1.611924550017006, -3.611254871967348};

    // max absolute value of the pressure difference
    const static double MAX_ABS_VAL_PRESSURE_DIFF = 106.0;
    const static double MAX_DISPLACEMENT_SYSTOLE  = 0.5; // papillary tips move this far at their peak
    
    // We require that the structures are associated with the finest level of
    // the patch hierarchy.
    const int level_num = hierarchy->getFinestLevelNumber();

    // Get the Lagrangian mesh.
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_num);
    const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();
    std::vector<LNode*> nodes = local_nodes;
    nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());

    // index without periodicity in Fourier series
    unsigned int k = (unsigned int) floor(current_time / (fourier_body_force->dt));
    
    // take periodic reduction                         
    unsigned int idx = k % (fourier_body_force->N_times);

    // current prescribed pressure difference
    const double pressure_mmHg = fourier_body_force->values[idx];

    // move compared to the current pressure difference
    // if the pressure is negative (higher ventricular pressure towards closure)
    const double displacement  = (pressure_mmHg < 0.0) ? MAX_DISPLACEMENT_SYSTOLE * abs(pressure_mmHg / MAX_ABS_VAL_PRESSURE_DIFF) : 0.0;


    // Loop over all Lagrangian mesh nodes and update the target point
    // positions.
    for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it){
        const LNode* const node = *it;
        IBTargetPointForceSpec* const force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec){
            // Here we update the position of the target point.
            //
            // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the total number of Lagrangian points)
            //        X_target     is the target position of the target point
            //        X_target(0)  is the x component of the target position
            //        X_target(1)  is the y component of the target position
            
            Point& X_target = force_spec->getTargetPointPosition();
            const int lag_idx = node->getLagrangianIndex();
            
            if (lag_idx == 0){
                // z coord down
                X_target(2) = LEFT_PAPILLARY[2]  - displacement;

            }
            else if (lag_idx == 1){
                // z coord down
                X_target(2) = RIGHT_PAPILLARY[2] - displacement;
            }
        }
    }

    return;
}// update_target_point_positions



inline double spring_function_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    /* 
    Compute force for collagen springs
    
    Zero force under compression
    
    Exponential growth until full recruitment,
    experimentally determined strain at which microscopic fibers align
     
    Linear after full recruitment
    
    All collagen assumed to have same modulus
    General parameters params[0] are normalized,
    but if a stiffer spring is requested then it is used here.
    
    params[0] must include the width element ds for converting forces to areas 
    */
    
    static const double a                    = 4643.4;       // Coeff of exponential term
    static const double b                    = 49.9643;      // Exponential rate
    static const double full_recruitment     = 0.145;             // Linear at strains larger than this
    static const double eta_collagen         = 32.5 * MPa_TO_CGS; // Linear slope, in barye = dynes/cm^2 = g cm/(s cm^2)
    static const double collagen_x_intercept = 0.125;             // Linear collagen part intercepts x axis at this strain
    static const double collagen_y_intercept = -collagen_x_intercept * eta_collagen; // Linear collagen part intercepts y axis at this stress
    static const double thickness            = 0.1; // cm, Thickness of leaflet tissue, for converting strains to forces
    
    // static const double ds, width element of the spring, also for converting strains to forces
    // included in params[0]
    
    const double kappa = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    

    
    // Compute the force
    if (E > full_recruitment){
        if ((lag_mastr_idx % 1000) == 0){
            std::cout << "On index " << lag_mastr_idx << "linear part, E = " << E << "\tindices = " << lag_mastr_idx << ", " << lag_slave_idx
            << ".\tEffective slope = " << kappa * thickness * eta_collagen
            << "\trest len = " << rest_len <<"\n";
        }
        return kappa * thickness * (eta_collagen*E + collagen_y_intercept);
    }
    else if (E > 0.0){
        if ((lag_mastr_idx % 5000) == 0){
            std::cout << "On index " << lag_mastr_idx << ",\t strain E = " << E << ",\texponential part, force = "
            << kappa * thickness * a * (exp(b*E) - 1) << "\tcoeff (using taylor series on E) = " << kappa * thickness * a * b
            << "\trest len = " << rest_len <<"\n";
        }
        return kappa * thickness * a * (exp(b*E) - 1);
    }
    else{
        return 0.0;
    }
    return 0.0;
} // spring_function_collagen

inline double deriv_spring_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx){
    // not implemented

    SAMRAI_MPI::abort();
    return 0.0;
} // deriv_spring_collagen





