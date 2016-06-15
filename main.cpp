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
#include <ibamr/IBInstrumentPanel.h>
#include <ibamr/IBStandardSourceGen.h>
#include <vector>
#include <queue>
#include <timing.h>
#include <boundary_condition_util.h>
#include <FeedbackForcer.h>
#include <FourierBodyForce.h>

// #define IMPLICIT_SOLVER
#ifdef IMPLICIT_SOLVER
     #include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#endif


#if defined(IBAMR_HAVE_SILO)
#include <silo.h>
#endif

void init_source_variables(vector<double>& Q_src, const vector<double>& P_src); 
void set_source_variables(vector<double>& Q_src, const vector<double>& P_src, double current_time, double dt, bool sink_on);
void update_rest_lengths(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* const l_data_manager, const double alpha);
                                 


//#define ENABLE_SOURCES 
//#define COMPARE_TO_ZERO_FLOW
#define DEBUG_OUTPUT 0 
#define ENABLE_INSTRUMENTS
#define FOURIER_SERIES_BC
#define DYNAMIC_BOUNDARY_STAB

// #define UPDATE_REST_LEN 

#define MMHG_TO_CGS 1333.22368
#define CGS_TO_MMHG 0.000750061683

#define SOURCE_ON_TIME       0.1
#define CONST_SRC_TIME       0.3
#define CONST_SRC_STRENGTH  93.0

#define MAX_STEP_FOR_CHANGE 1000

#define PI_DEFINED 3.1415926535897932384626433832795028841971693993751

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
        bool from_restart = false;
        if (argc >= 4)
            from_restart = true;

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
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);
        
        
        
        #ifdef ENABLE_SOURCES
            // set up the source
            Pointer<IBStandardSourceGen> ib_source; 
            std::vector<double> Q_src;
            std::vector<double> P_src;
            bool sink_on = false; 
        
            ib_source = new IBStandardSourceGen();
        
            if (ib_method_ops->hasFluidSources())
                pout << "ib_method_ops->hasFluidSources returned true before calling register\n" ;
            else 
                pout << "ib_method_ops->hasFluidSources returned FALSE before calling register\n" ;
        
            ib_method_ops->registerIBLagrangianSourceFunction(ib_source);

            if (ib_method_ops->hasFluidSources())
                pout << "ib_method_ops->hasFluidSources returned true AFTER calling register\n" ;
            else 
                pout << "ib_method_ops->hasFluidSources returned FALSE AFTER calling register\n" ;
                
        #endif
        
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
            
            // this fails here, need to get dt from input db or something 
            // dt = time_integrator->getMaximumTimeStepSize(); 
            
            dt = input_db->getDouble("DT"); 
            
            pout << "to constructor\n"; 
            fourier_series_data *fourier = new fourier_series_data("fourier_coeffs.txt", dt);  
            pout << "series data successfully built\n"; 
            // fourier->print_values(); 
        #endif
        
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
                VelocityBcCoefs *z_bdry_coeffs = new VelocityBcCoefs(fourier);
            
                u_bc_coefs[2] = z_bdry_coeffs;
            
                #ifdef DYNAMIC_BOUNDARY_STAB
            
                    // only for staggered grids 
                    if (solver_type == "STAGGERED"){
                        // always the Z component
                        Pointer<FeedbackForcer> bdry_dynamic_stab = new FeedbackForcer(z_bdry_coeffs, navier_stokes_integrator, patch_hierarchy);
                        time_integrator->registerBodyForceFunction(bdry_dynamic_stab);
                    }
                #endif
            
            #else 
                pout << "Using b.c. from file, no Fouier series created.\n"; 
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
        
        #ifdef ENABLE_SOURCES
            
            // more source stuff
            // can't put this before initializePatchHierarchy
            pout << "to initializeLevelData\n" ; 
            ib_source->initializeLevelData(patch_hierarchy, 0, 0.0, true, l_data_manager);
            std::cout << "cleared initializeLevelData" << std::endl;

            pout << "after initializePatchHierarchy and initializeLevelData\n" ; 
            
            // zero then set the source variables 
            // double dt_temp = 0.0;
            init_source_variables(ib_source->getSourceStrengths(finest_hier_level), ib_source->getSourcePressures(finest_hier_level));
            
            Q_src                 = ib_source->getSourceStrengths(finest_hier_level);
            const int num_sources = Q_src.size(); 
            const int source_idx  = 0; 
            const int sink_idx    = Q_src.size()-1; 
            
            // Stream to write-out pressure norm data 
            std::ofstream source_output_stream;
            if (!from_restart){
                if (SAMRAI_MPI::getRank() == 0){
                    source_output_stream.open("source_pressure_plot.m", ios_base::out | ios_base::trunc);
                    source_output_stream << "data = [" ; 
                }
            }
            // if we are restarting, we want to append to the file 
            // we are assuming that the pressure plot write was not interuppted here in the previous run 
            else{
                if (SAMRAI_MPI::getRank() == 0){
                    source_output_stream.open("source_pressure_plot.m", ios_base::out | ios_base::app);
                    source_output_stream << "data = [" ; 
                }
            }
            
        #endif


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


        // For debug, compare to zero flow 
        // Use an ifdef to avoid having things go out of scope 
        #ifdef COMPARE_TO_ZERO_FLOW
            
            // Create Eulerian initial condition specification objects.  These
            // objects also are used to specify exact solution values for error
            // analysis.
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
            
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction("p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
            
            // Setup data used to determine the accuracy of the computed solution.
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            
            const Pointer<Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
            const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();
            
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
            const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
            
            const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
            const Pointer<VariableContext> p_ctx = navier_stokes_integrator->getCurrentContext();
            
            const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
            const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
            //visit_data_writer->registerPlotQuantity("P error", "SCALAR", p_cloned_idx);
            
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx);
                patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx);
            }
            
            // Stream to write-out pressure norm data 
            if (!from_restart){
                std::ofstream pressure_plot;
                if (SAMRAI_MPI::getRank() == 0){
                    pressure_plot.open("pressure_plot.m", ios_base::out | ios_base::trunc);
                    pressure_plot << "pressure = [" ; 
                }
            }
            else{
                if (SAMRAI_MPI::getRank() == 0){
                    pressure_plot.open("pressure_plot.m", ios_base::out | ios_base::app);
                    pressure_plot << "pressure = [" ; 
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


        #ifdef UPDATE_REST_LEN
            
            // final rest length multiplier 
            double beta = 0.8; 
            
            // update over 200 times, 10 steps apart
            // for dt = 2.5e-5, update for 
            unsigned int num_spring_updates = 200; 
            unsigned int update_frequency   = 10;
            
            double alpha = pow(beta, 1.0 / ((double) num_spring_updates) ); 
            
            unsigned int updates_performed = 0; 
            
        #endif 


        #ifdef ENABLE_SOURCES
            if (from_restart){
                
                // do not update rest lengths if is from restart
                // if updates have not been completed, this is an error 
                if ( num_spring_updates * update_frequency > ((unsigned int) iteration_num) ){
                    pout << "Must have finished updating the spring rest lengths before the restart\n" ; 
                    SAMRAI_MPI::abort();  
                }
                
                // open the file, read and set
                char file_name[100]; 
                
                // this is already set above
                // iteration_num = time_integrator->getIntegratorStep();
                sprintf(file_name, "restart_src_vals_%d.txt", iteration_num); 
                
                FILE *file = fopen(file_name, "r");  
                if (!file){
                    pout << "problems opening restart file\n" ; 
                    SAMRAI_MPI::abort();
                }
                
                double source_restart_val; 
                double   sink_restart_val; 
                fscanf(file, "%lf %lf", &source_restart_val, &sink_restart_val); 
                
                Q_src[0]        = source_restart_val; 
                
                int sink_idx    = Q_src.size() - 1; 
                Q_src[sink_idx] = sink_restart_val;
                
                fclose(file); 
            }
        #endif 


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

            #ifdef UPDATE_REST_LEN
                if (!from_restart){
                    unsigned int it = (unsigned int) iteration_num; 
                    if( (updates_performed < num_spring_updates) && ((it % update_frequency) == 0)){
                    
                        update_rest_lengths(patch_hierarchy, l_data_manager, alpha); 
                        updates_performed++; 
                    } 
                }
            #endif 

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
            
            
            
            #ifdef ENABLE_SOURCES
                // reset source and sink for current time
                
                if ( loop_time >= SOURCE_ON_TIME ){ 
                    sink_on = true; 
                }
                
                set_source_variables(ib_source->getSourceStrengths(finest_hier_level), ib_source->getSourcePressures(finest_hier_level), loop_time, dt, sink_on);
                Q_src   = ib_source->getSourceStrengths(finest_hier_level);
                P_src   = ib_source->getSourcePressures(finest_hier_level);
                pout << "Source, Q = " << Q_src[source_idx] << ", P = " << P_src[source_idx] << "\n" ;
                pout << "Sink, Q = " << Q_src[sink_idx] << ", P = " << P_src[sink_idx] << "\n" ;
                            
            #endif

            
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
            
            #ifdef COMPARE_TO_ZERO_FLOW
                // Compute velocity and pressure error norms.
                const int coarsest_ln = 0;
                const int finest_ln = patch_hierarchy->getFinestLevelNumber();
                for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
                {
                    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                    if (!level->checkAllocated(u_cloned_idx)) level->allocatePatchData(u_cloned_idx);
                    if (!level->checkAllocated(p_cloned_idx)) level->allocatePatchData(p_cloned_idx);
                }
                
                pout << "\n" << "Computing error norms.\n\n";
                
                u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
                p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);
                
                HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
                hier_math_ops.setPatchHierarchy(patch_hierarchy);
                hier_math_ops.resetLevels(coarsest_ln, finest_ln);
                const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
                const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
                 
                // confirm that we have staggared grid data before outputting 
                Pointer<SideVariable<NDIM, double> > u_sc_var = u_var;
                if (u_sc_var)
                {
                    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
                    hier_sc_data_ops.subtract(u_cloned_idx, u_idx, u_cloned_idx);
                    pout << "Error in u at time " << loop_time << ":\n"
                    << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_cloned_idx, wgt_sc_idx) << "\n"
                    << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_cloned_idx, wgt_sc_idx) << "\n"
                    << "  max-norm: " << hier_sc_data_ops.maxNorm(u_cloned_idx, wgt_sc_idx) << "\n\n";
                }
            
                HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
                hier_cc_data_ops.subtract(p_cloned_idx, p_idx, p_cloned_idx);
                double pressure_norm = hier_cc_data_ops.maxNorm(p_cloned_idx, wgt_cc_idx);
                pout << "Norms of p at time " << loop_time - 0.5 * dt << ":\n" 
                //<< "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_cloned_idx, wgt_cc_idx) << "\n"
                //<< "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_cloned_idx, wgt_cc_idx) << "\n"
                << "  max-norm: " << pressure_norm * CGS_TO_MMHG << " mm Hg,    " << pressure_norm << " g/(cm s^2) \n\n";
                
                // write to file also 
                if (SAMRAI_MPI::getRank() == 0){
                    pressure_plot << pressure_norm << ", " ;
                    pressure_plot.flush(); 
                } 
                
                // see IB ex0 for more stuff here 
            #endif
            
            
            
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
                
                #ifdef ENABLE_SOURCES

                    // write source data for indices actually used
                    // write pressure data for all indices 
                    if (SAMRAI_MPI::getRank() == 0){
                    
                        source_output_stream << iteration_num << ", " << loop_time << ", " ;
                    
                        source_output_stream << Q_src[source_idx] << ", " << Q_src[sink_idx] ;
                    
                        for(int i=0; i<num_sources; i++){
                            source_output_stream << ", " << P_src[i]  ;
                        }
                    
                        source_output_stream << ";\n";                     
                        source_output_stream.flush(); 
                    }
            
                #endif
                
                #ifdef ENABLE_INSTRUMENTS
                    
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
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                
                #ifdef ENABLE_SOURCES
                    // write source and sink data from the current time 
                    std::ofstream source_data;
                    if (SAMRAI_MPI::getRank() == 0){
                        char file_name[100]; 
                        sprintf(file_name, "restart_src_vals_%d.txt", iteration_num); 
                        source_data.open(file_name, ios_base::out | ios_base::trunc);
                        source_data << std::setprecision(20) ;
                        source_data << Q_src[source_idx] << " " << Q_src[sink_idx] << "\n";  
                        source_data.close(); 
                    }
                #endif 
                
                
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
        
        
        #ifdef COMPARE_TO_ZERO_FLOW
            // Close the logging streams.
            if (SAMRAI_MPI::getRank() == 0){
                pressure_plot << "]; \n\n"; 
                pressure_plot.close();
            }
        #endif
        
        #ifdef ENABLE_SOURCES
            if (SAMRAI_MPI::getRank() == 0){
                source_output_stream << "]; \n\n"; 
                source_output_stream.close();
            }
        #endif
        
        #ifdef ENABLE_INSTRUMENTS
            if (SAMRAI_MPI::getRank() == 0){
                flux_output_stream << "]; \n\n"; 
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


void init_source_variables(vector<double>& Q_src, const vector<double>& P_src){
    /*
     Initializes source values Q to zero      
     */
    
    
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(Q_src.size() == P_src.size());
    // TBOX_ASSERT(Q_src.size() == 2);
#endif
    
    for(unsigned int i=0; i<Q_src.size(); i++){
        Q_src[i] = 0.0;
    }
    
    return;
}

void set_source_variables(vector<double>& Q_src, const vector<double>& P_src, double current_time, double dt, bool sink_on){
    /*
    Sets values Q. 
    
    All are considered sources, a sign swap is added for a sink.

    First value set to strength.
    
    Second uses resistor in series with inductor type model
    
    Resistance is set to support average flow at 100 mmHg
        R = (100 mmHg) / (5.6 L/min)
        R = 1429  (1 / (cm s^2)) in CGS units
    
    Attached to 5 mmHg pressure reservoir
    Set to near zero as crude model for venus pressure 
        P_0_sink = (80 mmHg) = 1333.22368 * 5 Barye
    
    Inertance (equivalent to inductance)
    Time scale is set assuming that dt = 2.5 * 10^(-5)
    L should be a few orders of magnitude larger than the
        L = 5
    
    Model:
        L (dQ/dt) + R Q = P - P_0_sink
    
    Discretized with implicit trapezoial rule in Q.
    Since P is computed at midpoints only, it is automatically time centered.
    (unclear whether this is actually correct, this is a step back in P)
    
        (L/dt)*(Q_n+1 - Q_n) + 0.5*(Q_n+1 + Q_n) = P - P_0_sink;
    
    This gives the NEGATIVE of the code value.
    */
    
    
//     const int num_src = Q_src.size();
    
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(Q_src.size() == P_src.size());
        // TBOX_ASSERT(Q_src.size() == 2);
#endif

    if (current_time < SOURCE_ON_TIME){
        Q_src[0] = 0.0;
    }    
    else if (current_time < CONST_SRC_TIME){
        Q_src[0] = CONST_SRC_STRENGTH * 0.5 * (1.0  - cos(PI_DEFINED * (current_time - SOURCE_ON_TIME) / (CONST_SRC_TIME - SOURCE_ON_TIME)));
    }
    else{
        Q_src[0] = CONST_SRC_STRENGTH; 
    }
    
    
    const static double R = 1429.0;
    const static double L = 5.0;
    const static double P_0_sink = MMHG_TO_CGS * 5.0;   
    
    
    if (sink_on){
    
        int sink_idx = Q_src.size() - 1; 
    
        double Q_sink_old = -Q_src[sink_idx];
        double Q_sink;
   
        Q_sink = ((L - 0.5*dt*R) * Q_sink_old + dt * (P_src[sink_idx] - P_0_sink)) / (L + 0.5*dt*R);
   
        Q_src[sink_idx] = -Q_sink;
    }

    return;
}







void update_rest_lengths(Pointer<PatchHierarchy<NDIM> > hierarchy, LDataManager* const l_data_manager, const double alpha){
    // Multiplies rest lengths of all local springs by multiplier
    
    
    // We require that the structures are associated with the finest level of
    // the patch hierarchy.
    const int level_num = hierarchy->getFinestLevelNumber();


    // Get the Lagrangian mesh.
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_num);
    
    const std::vector<LNode*>& nodes = l_mesh->getLocalNodes();
    
    //const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    //const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();
    // std::vector<LNode*> nodes = local_nodes;
    // nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());

    // Loop over all Lagrangian mesh nodes and update spring rest lengths 
    for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it){
    
        const LNode* const node = *it;
        IBSpringForceSpec*  spring_spec =  node->getNodeDataItem<IBSpringForceSpec>();
        
        if (spring_spec){
            vector<double>& params = spring_spec->getParameters()[0];
            params[1] *= alpha;  
        }
    }
    
    return;
}






