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

// Modified 2019, Alexander D. Kaiser

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
#include <CirculationModel_with_lv.h>
#include <FeedbackForcer.h>
#include <FourierBodyForce.h>
#include <cstdio>

#if defined(IBAMR_HAVE_SILO)
#include <silo.h>
#endif


typedef struct{

    // number of target points that move  
    int N_targets;      
    
    // vertex number for each target 
    int *vertex_idx;  
    
    // base coordinates 
    double *x_systole; 
    double *y_systole;
    double *z_systole;
    
    // maximum increment 
    double x_increment_systole_to_diastole; 
    double y_increment_systole_to_diastole; 
    double z_increment_systole_to_diastole; 

    // papillary movement follows these times 
    double t_diastole_start; 
    double t_diastole_full; 
    double t_systole_start; 
    double t_systole_full; 
    double t_cycle_length; 
    
    // papillary target point velocity
    // constant across all papillary points 
    double u_papillary[3]; 
    
} papillary_info; 

typedef struct{

    // number of target points that move  
    int N_vertices;      
    int N_times; 
    double dt_registration; 

    double t_smoothing; 

    // stored as 
    // positions[i + vertex_idx*3 + time_num*(3*N_vertices)]; 
    double *position; 
    
} prescribed_motion_info; 


inline double spring_function_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);
inline double deriv_spring_collagen(double R, const double* params, int lag_mastr_idx, int lag_slave_idx);

papillary_info* initialize_moving_papillary_info(string structure_name, 
                                                 fourier_series_data *fourier_series, 
                                                 LDataManager *l_data_manager); 

prescribed_motion_info* initialize_prescribed_motion_info(string structure_name, double t_cycle_length, double t_smoothing); 

void delete_prescribed_motion(prescribed_motion_info *motion_info); 

prescribed_motion_info* combine_prescribe_motion(prescribed_motion_info *motion_1, prescribed_motion_info *motion_2); 


void update_prescribed_motion_positions(Pointer<PatchHierarchy<NDIM> > hierarchy, 
                                        LDataManager* const l_data_manager, 
                                        const double current_time, 
                                        const double dt, 
                                        const double t_cycle_length, 
                                        fourier_series_data *fourier_series, 
                                        prescribed_motion_info* motion_info); 

void get_linear_interp_position(const double current_time, 
                                const double t_cycle_length, 
                                prescribed_motion_info* motion_info,
                                const int lag_idx,
                                double *position);

void get_smoothed_interp_position(const double current_time, 
                                  const double t_cycle_length, 
                                  prescribed_motion_info* motion_info,
                                  const int lag_idx,
                                  double *position); 

void print_prescribed_motion_summary(prescribed_motion_info *motion_info); 


#define DEBUG_OUTPUT 0 
#define ENABLE_INSTRUMENTS

// #define USE_CIRC_MODEL

#define MMHG_TO_CGS 1333.22368
#define CGS_TO_MMHG 0.000750061683
#define MPa_TO_CGS 1.0e7

#define SIM_WITH_LV

//#define DEBUG_PRESCRIBED_MOTION


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

        int n_restarts_written = 0;
        int max_restart_to_write = 15;
        if (input_db->keyExists("MAX_RESTART_TO_WRITE")){
            max_restart_to_write = input_db->getInteger("MAX_RESTART_TO_WRITE");
        }

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
        
            flux_output_stream.precision(14);
            flux_output_stream.setf(ios_base::scientific);
    
        
            // set below in loop
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
        // Fourier series first below 
        // vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));
                
        // This is needed to pull variables later
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        pout << "To Fourier series creation\n";
    
        dt = input_db->getDouble("DT");
        bool restart_circ_model = true;
        double start = input_db->getDouble("START_TIME");

        const bool use_circ_model = true; 
        // CirculationModel *circ_model   = new CirculationModel("circ_model", P_PA_0, P_LA_0, start, restart_circ_model);
        pout << "To constructor\n";

        // prescribing LV values at aorta outlet for now 
        fourier_series_data *fourier_aorta      = new fourier_series_data("fourier_coeffs_aorta.txt", dt);
        fourier_series_data *fourier_atrium     = new fourier_series_data("fourier_coeffs_atrium.txt", dt);
        fourier_series_data *fourier_ventricle  = new fourier_series_data("fourier_coeffs_ventricle.txt", dt);
        pout << "Series data successfully built\n";
    
        bool debug_series_values = false; 
        if (debug_series_values){
            pout << "fourier_aorta = [\n"; 
            fourier_aorta->print_values(); 
            pout << "];\n"; 
            pout << "fourier_atrium = [\n"; 
            fourier_atrium->print_values(); 
            pout << "];\n";
            pout << "fourier_ventricle = [\n"; 
            fourier_ventricle->print_values(); 
            pout << "];\n";
        }


        // scaled cycle length for this patient 
        double t_cycle_length = input_db->getDouble("CYCLE_DURATION");

        // start at systole, which is this far into the Fourier series 
        double t_offset_start_bcs_unscaled = input_db->getDouble("T_OFFSET_START_BCS_UNSCALED"); // starts this 

        // start in physical time with relation to Fourier series 
        double t_offeset_start = t_offset_start_bcs_unscaled * (t_cycle_length / fourier_atrium->L); 

        std::string structure_name = input_db->getString("NAME"); 
                
        // this just sets necessary data structures to zero in this code 
        papillary_info* papillary = initialize_moving_papillary_info(structure_name, fourier_atrium, l_data_manager); 

        std::string structure_name_LV = input_db->getString("NAME_LV"); 

        pout << "structure_name_LV = " << structure_name_LV << "\n"; 

        // set up the ventricle and valve skeleton motion 
        double t_smoothing = 4.0e-2; 
        prescribed_motion_info* motion_info_LV     = initialize_prescribed_motion_info(structure_name_LV, t_cycle_length, t_smoothing); 
        prescribed_motion_info* motion_info_mitral = initialize_prescribed_motion_info(structure_name, t_cycle_length, t_smoothing); 
        prescribed_motion_info* motion_info        = combine_prescribe_motion(motion_info_LV, motion_info_mitral); 
        delete_prescribed_motion(motion_info_LV); 
        delete_prescribed_motion(motion_info_mitral);         


        string aorta_vertices_file_name  = "aorta_bdry.vertex"; 
        string atrium_vertices_file_name = "atrium_bdry.vertex"; 

        CirculationModel_with_lv * circ_model_with_lv = new CirculationModel_with_lv(fourier_aorta, 
                                                                                     fourier_atrium, 
                                                                                     fourier_ventricle,
                                                                                     aorta_vertices_file_name,
                                                                                     atrium_vertices_file_name,
                                                                                     t_cycle_length,
                                                                                     t_offset_start_bcs_unscaled,
                                                                                     time_integrator->getIntegratorTime()); // current time here 

        // Create Eulerian boundary condition specification objects.
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (int d = 0; d < NDIM; ++d){
            u_bc_coefs[d] = new VelocityBcCoefs_lv_aorta(d, circ_model_with_lv);
        }
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        // flow straightener at boundary 
        Pointer<FeedbackForcer> feedback_forcer = new FeedbackForcer(navier_stokes_integrator, patch_hierarchy, circ_model_with_lv);
        time_integrator->registerBodyForceFunction(feedback_forcer);

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
                
            #endif
            
        }
 
        
        double dt_original = time_integrator->getMaximumTimeStepSize();


        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        dt = 0.0;
        double dt_prev = 0.0;
        bool prev_step_initialized = false;
        
        get_timestamp(&time2_total);
        double total_init_time = timestamp_diff_in_seconds(time1_total, time2_total);
        pout << "\n\nTotal initialization time = " << total_init_time << "\n\n"; 
        
        #ifdef DEBUG_PRESCRIBED_MOTION
            dt = 5.0e-4; 
            double max_t = 2*t_cycle_length;
            for(double current_time=0.0; current_time<max_t; current_time+=dt){
                update_prescribed_motion_positions(patch_hierarchy, 
                                                  l_data_manager, 
                                                  current_time, 
                                                  dt, 
                                                  t_cycle_length, 
                                                  fourier_atrium, 
                                                  motion_info);
            }

            // finalize 
            dt = 0.0; 
            update_prescribed_motion_positions(patch_hierarchy, 
                                  l_data_manager, 
                                  t_cycle_length, 
                                  dt, 
                                  t_cycle_length, 
                                  fourier_atrium, 
                                  motion_info);

            pout << "Abort after motion debug\n"; 
            SAMRAI_MPI::abort();

        #endif



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
            
            // update circ model  
            {
                Pointer<hier::Variable<NDIM> > U_var = navier_stokes_integrator->getVelocityVariable();
                Pointer<hier::Variable<NDIM> > P_var = navier_stokes_integrator->getPressureVariable();
                Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
                const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, current_ctx);
                const int P_current_idx = var_db->mapVariableAndContextToIndex(P_var, current_ctx);
                Pointer<HierarchyMathOps> hier_math_ops = navier_stokes_integrator->getHierarchyMathOps();
                const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
                const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
                circ_model_with_lv->advanceTimeDependentData(dt, patch_hierarchy, U_current_idx, P_current_idx, wgt_cc_idx, wgt_sc_idx);
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

             
            // LV motion 
            update_prescribed_motion_positions(patch_hierarchy, 
                                              l_data_manager, 
                                              loop_time, 
                                              dt, 
                                              t_cycle_length, 
                                              fourier_atrium, 
                                              motion_info); 
            

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
            
                circ_model_with_lv->set_Q_mitral(-flux_valve_ring[0]); 

                // body_force->d_flux_z = flux_valve_ring[0]; 
            
            #endif
                        
            
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                
                n_restarts_written++;
                
                if (n_restarts_written > max_restart_to_write){
                    if (SAMRAI_MPI::getRank() == 0){
                        std::ofstream controlled_stop_stream;
                        controlled_stop_stream.open("controlled_stop.txt", ios_base::out | ios_base::trunc);
                        controlled_stop_stream << "stopped\n" ;
                        controlled_stop_stream.close();
                    }
                    SAMRAI_MPI::barrier();
                    SAMRAI_MPI::abort();
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
                
                flux_output_stream << "fig = figure;\n plot(data(:,1), -data(:,2), 'k');\n";
                flux_output_stream << "hold on;\n";
                flux_output_stream << "dt = " << dt_original << "; \n"; 
                flux_output_stream << "net_flux = dt*cumsum(-data(:,2));\n "; 
                flux_output_stream << "plot(data(:,1), net_flux, '--k');\n";
                flux_output_stream << "xlabel('t');\n ylabel('ml/s, ml');\n";
                flux_output_stream << "legend('flux', 'net flux', 'Location', 'NorthWest')\n";
                flux_output_stream << "plot(data(:,1), 0*data(:,2), ':k');\n";
                flux_output_stream << "printfig(fig,'flux.eps');\n"; 
                
                flux_output_stream.close();
            }
        #endif
        
        #ifdef USE_CIRC_MODEL
            circ_model->write_plot_code(); 
        #endif 
        
        if (SAMRAI_MPI::getRank() == 0){
            std::ofstream done_stream;
            done_stream.open("done.txt", ios_base::out | ios_base::trunc);
            done_stream << "done\n" ;
            done_stream.close(); 
        }

        delete_prescribed_motion(motion_info); 

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main



papillary_info* initialize_moving_papillary_info(string structure_name, 
                                                 fourier_series_data *fourier_series, 
                                                 LDataManager *l_data_manager){
    // reads file structure_name.papillary 
    // and initializes information for moving papillary 
    
    papillary_info* papillary = new papillary_info;
    
    string full_name = structure_name + string(".papillary"); 
    
    std::cout << "name = " << full_name << "\n"; 
    
    ifstream papillary_file(full_name.c_str());
    
    // first is number of points 
    papillary_file >> papillary->N_targets; 
    
    papillary->vertex_idx = new int[papillary->N_targets];  
    papillary->x_systole  = new double[papillary->N_targets]; 
    papillary->y_systole  = new double[papillary->N_targets]; 
    papillary->z_systole  = new double[papillary->N_targets]; 
    
    // next gets x,y,z increments 
    papillary_file >> papillary->x_increment_systole_to_diastole; 
    papillary_file >> papillary->y_increment_systole_to_diastole; 
    papillary_file >> papillary->z_increment_systole_to_diastole;
    
    papillary_file >> papillary->t_diastole_start; 
    papillary_file >> papillary->t_diastole_full; 
    papillary_file >> papillary->t_systole_start; 
    papillary_file >> papillary->t_systole_full; 
    papillary_file >> papillary->t_cycle_length;
    
    if (papillary->t_diastole_start != 0.0){
        std::cout << "Diastole assumed to start at zero for current movement.\n"; 
        SAMRAI_MPI::abort();
    }
    
    std::cout << "N_targets = " << papillary->N_targets << "\n"; 
    std::cout << "increment = " << papillary->x_increment_systole_to_diastole << " " << papillary->y_increment_systole_to_diastole << " " << papillary->z_increment_systole_to_diastole << "\n"; 
    
    for (int i=0; i<papillary->N_targets; i++){
        papillary_file >> papillary->vertex_idx[i];  
        papillary_file >> papillary->x_systole[i]; 
        papillary_file >> papillary->y_systole[i]; 
        papillary_file >> papillary->z_systole[i];
        
        std::cout << "idx, coords = " << papillary->vertex_idx[i]  << " " <<  papillary->x_systole[i]  << " " << papillary->y_systole[i]  << " " << papillary->z_systole[i] << "\n";
    }
    
    papillary->u_papillary[0] = 0.0; 
    papillary->u_papillary[1] = 0.0; 
    papillary->u_papillary[2] = 0.0; 
    
    l_data_manager->initialize_movement_info(papillary->N_targets, papillary->vertex_idx, papillary->u_papillary); 
        
    std::cout << "Movement information in initialize\n"; 
    l_data_manager->print_movement_info(); 
        
    return papillary; 
}  

prescribed_motion_info* initialize_prescribed_motion_info(string structure_name, double t_cycle_length, double t_smoothing){

    prescribed_motion_info* motion_info = new prescribed_motion_info; 

    motion_info->t_smoothing = t_smoothing; 

    motion_info->N_times = 0; 

    // find number of files, number of vertices per file and total 
    for(int file_num=0; ; file_num++){
        char buffer[100];
        sprintf(buffer, "%02d", file_num);  
        string full_name = structure_name + string(buffer) + string(".vertex"); 

        pout << "looking for file " << full_name << "\n"; 

        ifstream current_file(full_name.c_str());        

        if (current_file){
            pout << "found file\n"; 
            (motion_info->N_times)++; 

            if (file_num == 0){
                current_file >> motion_info->N_vertices; 
            }
            else{
                int N_vertices_tmp; 
                current_file >> N_vertices_tmp; 
                if (N_vertices_tmp != motion_info->N_vertices){
                    pout << "Inconsistent number of vertices.\n"; 
                    SAMRAI_MPI::abort();
                }
            }
            current_file.close(); 

        }
        else{
            // done reading here 
            pout << "Did not find file, moving on\n"; 
            break; 
        }
    }

    motion_info->dt_registration = t_cycle_length / motion_info->N_times; 

    if (motion_info->dt_registration < (motion_info->t_smoothing)){
        pout << "Smoothing too large for existing registration time step\n"; 
        SAMRAI_MPI::abort();  
    } 

    int total_scalars = 3 * motion_info->N_vertices * motion_info->N_times; 
    motion_info->position = new double[total_scalars]; 


    for(int file_num=0; file_num < motion_info->N_times; file_num++){
        char buffer[100];
        sprintf(buffer, "%02d", file_num);  
        string full_name = structure_name + string(buffer) + string(".vertex"); 
        ifstream current_file(full_name.c_str());        

        if (current_file){

            int N_vertices_tmp; 
            current_file >> N_vertices_tmp; 
            if (N_vertices_tmp != motion_info->N_vertices){
                pout << "Inconsistent number of vertices.\n"; 
                SAMRAI_MPI::abort();
            }
                        
            for(int vertex_idx=0; vertex_idx<motion_info->N_vertices; vertex_idx++){
                for(int i=0; i<3; i++){
                    current_file >> motion_info->position[i + vertex_idx*3 + file_num*(3*motion_info->N_vertices)]; 
                }
            }
            current_file.close(); 
        }
        else{
            pout << "Could not open needed file\n"; 
            SAMRAI_MPI::abort(); 
        }

    }

    print_prescribed_motion_summary(motion_info); 

    return motion_info; 
}

void print_prescribed_motion_summary(prescribed_motion_info *motion_info){

    pout << "\n\n"; 
    pout << "N_vertices = " << motion_info->N_vertices << "\n";       
    pout << " N_times = " << motion_info->N_times << "\n"; 
    pout << " dt_registration = " << motion_info->dt_registration << "\n";  
    pout << "t_smoothing = " << motion_info->t_smoothing << "\n";  

    pout << "First few vertices at a few time steps:\n\n"; 

    for (int time_num=0; time_num<motion_info->N_times; time_num++){
        pout << "time_num = " << time_num << " time (s) = " << time_num * motion_info->dt_registration << "\n"; 

        for (int vertex_idx=0; vertex_idx<5; vertex_idx++){
            pout << "vertex_idx = " << vertex_idx << ", pos = ("; 

            for (int i=0; i<3; i++){
                pout << motion_info->position[i + vertex_idx*3 + time_num*(3*motion_info->N_vertices)] << " ";
            }

            pout << ")\n"; 
        }
    }
}

void delete_prescribed_motion(prescribed_motion_info *motion_info){
    delete motion_info->position; 
    delete motion_info; 
}

prescribed_motion_info* combine_prescribe_motion(prescribed_motion_info *motion_1, prescribed_motion_info *motion_2){
    // combine multiple motion info files into one 

    // some light error checking for compatible data structures here 

    if (motion_1->N_times != motion_2->N_times){
        pout << "Number of time steps not equal while combining struct\n"; 
        SAMRAI_MPI::abort(); 
    }
    if (motion_1->dt_registration != motion_2->dt_registration){
        pout << "registration time steps not equal while combining struct\n"; 
        SAMRAI_MPI::abort(); 
    }
    if (motion_1->t_smoothing != motion_2->t_smoothing){
        pout << "smoothing not equal while combining struct\n"; 
        SAMRAI_MPI::abort(); 
    }

    prescribed_motion_info* motion_combined = new prescribed_motion_info; 

    motion_combined->N_times         = motion_1->N_times; 
    motion_combined->dt_registration = motion_1->dt_registration; 
    motion_combined->t_smoothing     = motion_1->t_smoothing; 

    motion_combined->N_vertices      = motion_1->N_vertices + motion_2->N_vertices; 

    int total_scalars = 3 * motion_combined->N_vertices * motion_combined->N_times; 
    motion_combined->position = new double[total_scalars]; 

    for(int file_num=0; file_num < motion_1->N_times; file_num++){

        int int_idx_combined = 0; 
        for(int vertex_idx=0; vertex_idx<motion_1->N_vertices; int_idx_combined++, vertex_idx++){
            for(int i=0; i<3; i++){
                motion_combined->position[i + int_idx_combined*3 + file_num*(3*motion_combined->N_vertices)] = motion_1->position[i + vertex_idx*3 + file_num*(3*motion_1->N_vertices)]; 
            }
        }

        // no re-init of combined index 
        for(int vertex_idx=0; vertex_idx<motion_2->N_vertices; int_idx_combined++, vertex_idx++){
            for(int i=0; i<3; i++){
                motion_combined->position[i + int_idx_combined*3 + file_num*(3*motion_combined->N_vertices)] = motion_2->position[i + vertex_idx*3 + file_num*(3*motion_2->N_vertices)]; 
            }
        }

    }

    return motion_combined; 
}


void update_prescribed_motion_positions(Pointer<PatchHierarchy<NDIM> > hierarchy, 
                                        LDataManager* const l_data_manager, 
                                        const double current_time, 
                                        const double dt, 
                                        const double t_cycle_length, 
                                        fourier_series_data *fourier_series, 
                                        prescribed_motion_info* motion_info){


    // We require that the structures are associated with the finest level of
    // the patch hierarchy.
    const int level_num = hierarchy->getFinestLevelNumber();

    #ifdef DEBUG_PRESCRIBED_MOTION
        static std::ofstream plot_motion_debug; 
        static int idx_to_write = 0; 
        if (current_time == 0.0){
            plot_motion_debug.open("plot_motion_debug.m", ios_base::out | ios_base::trunc);
            plot_motion_debug.precision(14);
            plot_motion_debug << "data = [";
        }
    #endif 

    #ifdef DEBUG_PRESCRIBED_MOTION
        // pout << "checking dt... dt = " << dt << "\n"; 
        if (dt == 0.0){
            pout << "hit finalize in debug motion file output\n"; 
            plot_motion_debug << "]; figure;\n"; 
            plot_motion_debug << "plot(data(:,1), data(:,2), 'k');\n"; 
            plot_motion_debug << "figure;\n plot(data(:,1), data(:,3), 'k');\n";
            plot_motion_debug << "figure;\n plot(data(:,1), data(:,4), 'k');\n";
            plot_motion_debug << "figure;\n plot(data(:,1), data(:,5), 'k'); title('fraction to next step');\n";
            plot_motion_debug.close(); 
            return; 
        }
    #endif 


    // Get the Lagrangian mesh.
    Pointer<LMesh> l_mesh = l_data_manager->getLMesh(level_num);
    const std::vector<LNode*>& local_nodes = l_mesh->getLocalNodes();
    const std::vector<LNode*>& ghost_nodes = l_mesh->getGhostNodes();
    std::vector<LNode*> nodes = local_nodes;
    nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());

    double t_reduced = current_time - t_cycle_length * floor(current_time / t_cycle_length); 

    unsigned int min_step_motion = floor(t_reduced / motion_info->dt_registration); 
    unsigned int next_step_motion = (min_step_motion+1) % motion_info->N_times; 

    double min_step_time  = min_step_motion * motion_info->dt_registration; 
    double next_step_time = next_step_motion * motion_info->dt_registration; 

    double fraction_to_next_step = (t_reduced - min_step_time) / motion_info->dt_registration; 

    // // for smoothing via averaging 
    bool smoothing_on = false; 
    int smoothing_steps_per_side; 
    int smoothing_steps_total; 

    if (fabs(min_step_time - t_reduced) < motion_info->t_smoothing/2){
        // if we are right by the minimum time, then smoothing is on with the prevoius value         
        // do not smooth on initial step, only if we are above 
        if (current_time > (motion_info->t_smoothing)){
            smoothing_on = true; 
        }
    }
    else if (fabs(next_step_time - t_reduced) < motion_info->t_smoothing/2){
        // if we are right by the maximum time, then smoothing is on with the next, future value 
        smoothing_on = true; 
    }
    else if ((next_step_time == 0) && (fabs(next_step_time - (t_reduced - t_cycle_length)) < motion_info->t_smoothing/2)){
        // if we are right by the maximum time, counting a periodic wrap
        // this is when next_step_time has been reduced mod N and 
        smoothing_on = true; 
    }

    if (smoothing_on){
        smoothing_steps_per_side = floor( (motion_info->t_smoothing/2) / dt); 
        smoothing_steps_total    = smoothing_steps_per_side*2 + 1; 
    }

    bool print_summary = false; 
    int n_valve_printed = 0; 
    if (print_summary){
        pout << "motion update debug summary:\n"; 
        pout << "t_reduced = " << t_reduced << "\n"; 
        pout << "min_step_motion = " << min_step_motion << "\n"; 
        pout << "next_step_motion = " << next_step_motion << "\n"; 
        pout << "min_step_time = " << min_step_time << "\n"; 
        pout << "fraction_to_next_step = " << fraction_to_next_step << "\n"; 
    }

    if ((fraction_to_next_step < 0.0) || (fraction_to_next_step > 1.0)){
        pout << "fraction_to_next_step is not in allowable range\n";  
        SAMRAI_MPI::abort(); 
    }

    // Loop over all Lagrangian mesh nodes and update the target point
    // positions.
    for (std::vector<LNode*>::const_iterator it = nodes.begin(); it != nodes.end(); ++it){
        const LNode* const node = *it;
        IBTargetPointForceSpec* const force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
        
        if (force_spec && (force_spec->getStiffness() > 0.0)){ 
            // Here we update the position of the target point.
            //
            // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the total number of Lagrangian points)
            //        X_target     is the target position of the target point
            //        X_target(0)  is the x component of the target position
            //        X_target(1)  is the y component of the target position
            
            Point& X_target = force_spec->getTargetPointPosition();
            const int lag_idx = node->getLagrangianIndex();
            
            if (lag_idx < motion_info->N_vertices){

                double position[3] = {0.0, 0.0, 0.0};  

                // X_target(0) = (1.0 - fraction_to_next_step) * motion_info->position[0 + 3*lag_idx + min_step_motion *(3*motion_info->N_vertices)] + 
                //               (      fraction_to_next_step) * motion_info->position[0 + 3*lag_idx + next_step_motion*(3*motion_info->N_vertices)]; 
                // X_target(1) = (1.0 - fraction_to_next_step) * motion_info->position[1 + 3*lag_idx + min_step_motion *(3*motion_info->N_vertices)] + 
                //               (      fraction_to_next_step) * motion_info->position[1 + 3*lag_idx + next_step_motion*(3*motion_info->N_vertices)];
                // X_target(2) = (1.0 - fraction_to_next_step) * motion_info->position[2 + 3*lag_idx + min_step_motion *(3*motion_info->N_vertices)] + 
                //               (      fraction_to_next_step) * motion_info->position[2 + 3*lag_idx + next_step_motion*(3*motion_info->N_vertices)];

                if (smoothing_on){
                    get_smoothed_interp_position(current_time, 
                                                t_cycle_length, 
                                                motion_info,
                                                lag_idx,
                                                position);
                }
                else{
                    // take the current position from linear interpolation 
                    get_linear_interp_position(current_time, 
                                               t_cycle_length, 
                                               motion_info,
                                               lag_idx,
                                               position); 
                }

                for (int component=0; component<3; component++){
                    X_target(component) = position[component]; 
                } 

                if (print_summary){
                    if (lag_idx < 6) {
                        pout << "lag idx = " << lag_idx << ", pos = (" << X_target(0) << " " << X_target(1) << " " << X_target(2) << ")\n"; 
                    }
                    if (lag_idx == 78296){
                        pout << "lag idx = " << lag_idx << "(mitral idx 9624), pos = (" << X_target(0) << " " << X_target(1) << " " << X_target(2) << ")\n"; 
                    }
                }

                #ifdef DEBUG_PRESCRIBED_MOTION
                    if (lag_idx == idx_to_write){
                        plot_motion_debug << current_time << ", " << X_target(0) << ", " << X_target(1) << ", " << X_target(2) << ", " << fraction_to_next_step << "\n";  
                    }
                #endif 
            }

            // if(print_summary){
            //     if (lag_idx > motion_info->N_vertices){
            //         if(n_valve_printed < 6){
            //             pout << "lag idx = " << lag_idx << "mitral file idx = " << lag_idx - motion_info->N_vertices << ", pos = (" << X_target(0) << " " << X_target(1) << " " << X_target(2) << ")\n";
            //             n_valve_printed++; 
            //         } 
            //     }
            // }

        }
    }

}


void get_linear_interp_position(const double current_time, 
                                const double t_cycle_length, 
                                prescribed_motion_info* motion_info,
                                const int lag_idx,
                                double *position){

    double t_reduced = current_time - t_cycle_length * floor(current_time / t_cycle_length); 

    unsigned int min_step_motion = floor(t_reduced / motion_info->dt_registration); 
    unsigned int next_step_motion = (min_step_motion+1) % motion_info->N_times; 

    double min_step_time  = min_step_motion * motion_info->dt_registration; 

    double fraction_to_next_step = (t_reduced - min_step_time) / motion_info->dt_registration; 

    for(int i=0; i<3; i++){
        position[i] = (1.0 - fraction_to_next_step) * motion_info->position[i + 3*lag_idx + min_step_motion *(3*motion_info->N_vertices)] + 
                      (      fraction_to_next_step) * motion_info->position[i + 3*lag_idx + next_step_motion*(3*motion_info->N_vertices)]; 
    }

}



void get_smoothed_interp_position(const double current_time, 
                                  const double t_cycle_length, 
                                  prescribed_motion_info* motion_info,
                                  const int lag_idx,
                                  double *position){

    // takes convolution over the specified interval with a constant 

    // reduced mod cycles 
    double t_reduced = current_time - t_cycle_length * floor(current_time / t_cycle_length); 

    double dt_reg = motion_info->dt_registration; 

    // bounds of integral 
    double a = t_reduced - motion_info->t_smoothing/2.0; 
    double b = t_reduced + motion_info->t_smoothing/2.0; 

    // there is one time point that lies in the domain of integration 
    unsigned int center_step_motion = (unsigned int) floor(b / dt_reg);
    double t_cent = center_step_motion * dt_reg; 

    // idx gets periodic reduction
    center_step_motion %= motion_info->N_times; 

    // this time point is below the domain of integration
    unsigned int prev_step_motion; 
    double t_prev = t_cent - dt_reg; 
    // handle periodicity on indices 
    if (center_step_motion == 0){
        // index gets a periodic wrap 
        prev_step_motion   =  motion_info->N_times - 1;  

        // actual value is negative in time relative to current position 
        //t_prev = -dt_reg;
    }
    else{
        prev_step_motion = center_step_motion - 1; 
        //t_prev = prev_step_motion * dt_reg; 
    }

    // handle periodicity on indices 
    unsigned int next_step_motion = (center_step_motion+1) % motion_info->N_times; 
    // no periodicity on time 
    double t_next = t_cent + dt_reg; 

    // coefficients for each of the three time points 
    // integrate the piece-wise linear interpolant to get the following closed form on the weights 
    // double weight_prev   = t_cent*t_cent/2 - t_cent*a + a*a/2; 
    // double weight_center = t_cent*t_cent/2 - t_prev*t_cent - a*a/2 + t_prev*a + 
    //                        t_next*b -b*b/2 - t_next*t_cent + t_cent*t_cent/2; 
    // double weight_next   = b*b/2 - t_cent*b + t_cent*t_cent/2; 

    double weight_prev   =  pow(t_cent - a     ,2)/2; 
    double weight_center =  pow(t_cent - t_prev,2)/2 - pow(a      - t_prev,2)/2 + 
                           -pow(t_next - b     ,2)/2 + pow(t_next - t_cent,2)/2; 
    double weight_next   =  pow(b      - t_cent,2)/2; 

    // dt_reg comes from coefficients in interpolation function 
    // t_smoothing comes from that the integral is an average 
    weight_prev   /= (dt_reg * motion_info->t_smoothing); 
    weight_center /= (dt_reg * motion_info->t_smoothing); 
    weight_next   /= (dt_reg * motion_info->t_smoothing); 

    // if (lag_idx == 0){
    //     pout << "current_time = "<< current_time << " t_reduced = " << t_reduced << "\n"; 
    //     pout << "a,b = " << a << ", " << b << "\n"; 
    //     pout << "prev_step_motion = " << prev_step_motion << " t_prev_step = " << t_prev << "\n";
    //     pout << "center_step_motion = " << center_step_motion << " t_center = " << t_cent << "\n";
    //     pout << "next_step_motion = " << next_step_motion << " t_next_step = " << t_next << "\n";
    //     pout << "weights = " << weight_prev << ", " << weight_center << ", " << weight_next << "\n";
    //     pout << "\n";  
    // }


    for(int i=0; i<3; i++){
        position[i] = weight_prev   * motion_info->position[i + 3*lag_idx + prev_step_motion  *(3*motion_info->N_vertices)] + 
                      weight_center * motion_info->position[i + 3*lag_idx + center_step_motion*(3*motion_info->N_vertices)] +
                      weight_next   * motion_info->position[i + 3*lag_idx + next_step_motion  *(3*motion_info->N_vertices)];  
    }


}






void update_target_point_positions(Pointer<PatchHierarchy<NDIM> > hierarchy, 
                                   LDataManager* const l_data_manager, 
                                   const double current_time, 
                                   const double dt, 
                                   fourier_series_data *fourier_series, 
                                   papillary_info *papillary){
                                   
    // Requires to have pressure difference as a Fourier series
    // Atrial pressure is positive if higher
    // so positive pressure drives forward flow

    // quick return at beginning so that things do not move in discontinuous manner
    // magic number here, 0.1255 is the location in time of the max fwd pressure in diastole 
    //if (current_time == papillary->max_p_time)
    //if (current_time == 0.0)
    //    return;

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
/*    unsigned int k = (unsigned int) floor(current_time / (fourier_series->dt));
    
    // take periodic reduction                         
    unsigned int idx = k % (fourier_series->N_times);

    // current prescribed pressure difference
    const double pressure_mmHg = fourier_series->values[idx];
*/ 
    // move compared to the current pressure difference
    // if the pressure is negative (higher ventricular pressure towards closure)
    // double power = 1.0 / 10.0;
    double frac_to_diastole;
    double u_target[3]; 
    
    // papillary movement follows these times 
    double t_reduced = current_time - papillary->t_cycle_length * floor(current_time/ (papillary->t_cycle_length)); 
    // std::cout << "t = " << current_time << ", t_reduced = " << t_reduced << "\n"; 
    
    if (t_reduced < papillary->t_diastole_full){
        
        #ifdef C1_MOVEMENT
        
            frac_to_diastole = 0.5 * (1.0 - cos(M_PI * t_reduced/papillary->t_diastole_full));
            
            const double deriv_unscaled = M_PI/(2.0 * papillary->t_diastole_full) * sin(M_PI * t_reduced/papillary->t_diastole_full); 
        
            // constant velocity in linear movement 
            u_target[0] = papillary->x_increment_systole_to_diastole * deriv_unscaled;
            u_target[1] = papillary->y_increment_systole_to_diastole * deriv_unscaled;
            u_target[2] = papillary->z_increment_systole_to_diastole * deriv_unscaled;
            
        #else
        
            // linear movement 
            frac_to_diastole = t_reduced / papillary->t_diastole_full; 
        
            // constant velocity in linear movement 
            u_target[0] = papillary->x_increment_systole_to_diastole / papillary->t_diastole_full; 
            u_target[1] = papillary->y_increment_systole_to_diastole / papillary->t_diastole_full; 
            u_target[2] = papillary->z_increment_systole_to_diastole / papillary->t_diastole_full;  
                  
        #endif
        
    }
    else if  (t_reduced < papillary->t_systole_start){
        frac_to_diastole = 1.0; 
        
        // still in diastole 
        u_target[0] = 0.0;
        u_target[1] = 0.0;
        u_target[2] = 0.0;
        
    }
    else if  (t_reduced < papillary->t_systole_full){
    
        #ifdef C1_MOVEMENT
        
            frac_to_diastole = 0.5 * (cos(M_PI * (t_reduced - papillary->t_systole_start)/(papillary->t_systole_full - papillary->t_systole_start)) + 1.0); 
            
            const double deriv_unscaled = -0.5 *  sin(M_PI * (t_reduced - papillary->t_systole_start)/(papillary->t_systole_full - papillary->t_systole_start)) * (M_PI/(papillary->t_systole_full - papillary->t_systole_start)); 
        
            // constant velocity in linear movement 
            u_target[0] = papillary->x_increment_systole_to_diastole * deriv_unscaled;
            u_target[1] = papillary->y_increment_systole_to_diastole * deriv_unscaled;
            u_target[2] = papillary->z_increment_systole_to_diastole * deriv_unscaled;
        
        #else
        
            const double slope     = -1.0/(papillary->t_systole_full - papillary->t_systole_start); 
            const double intercept = papillary->t_systole_full/(papillary->t_systole_full - papillary->t_systole_start); 
            frac_to_diastole       = slope * t_reduced + intercept; 
        
            u_target[0] = -papillary->x_increment_systole_to_diastole / (papillary->t_systole_full - papillary->t_systole_start);
            u_target[1] = -papillary->y_increment_systole_to_diastole / (papillary->t_systole_full - papillary->t_systole_start);
            u_target[2] = -papillary->z_increment_systole_to_diastole / (papillary->t_systole_full - papillary->t_systole_start);
        
        #endif 
        
    }
    else{ 
        // full systole 
        frac_to_diastole = 0.0; 
                
        // still in full systole  
        u_target[0] = 0.0;
        u_target[1] = 0.0;
        u_target[2] = 0.0;
    }
    
    // update velocity with data manager 
    l_data_manager->set_movement_velocity(u_target); 
    
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
            
            // loop over struct, little wasteful but not that many 
            for (int i=0; i<(papillary->N_targets); i++){
                if (lag_idx == (papillary->vertex_idx[i])){
                    X_target(0) = papillary->x_systole[i] + frac_to_diastole * papillary->x_increment_systole_to_diastole;
                    X_target(1) = papillary->y_systole[i] + frac_to_diastole * papillary->y_increment_systole_to_diastole;
                    X_target(2) = papillary->z_systole[i] + frac_to_diastole * papillary->z_increment_systole_to_diastole;
                }            
            }
            
        }
    }

    // std::cout << "Movement information in update_target_point_positions:\n"; 
    // l_data_manager->print_movement_info(); 

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
    //static const double thickness            = 0.1; // cm, Thickness of leaflet tissue, for converting strains to forces
    
    // static const double ds, width element of the spring, also for converting strains to forces
    // included in params[0]
    
    const double kappa    = params[0];
    const double rest_len = params[1];
    
    // Strain, dimension 1
    const double E = R/rest_len - 1.0;
    
    // Compute the force
    if (E > full_recruitment){
        /*if ((lag_mastr_idx % 2500) == 0){
            std::cout << "Affine. (idx,nbr) = (" << lag_mastr_idx << ", " <<  lag_slave_idx
                      << "\tE = " << E
                      << "\tF = " << kappa * (eta_collagen*E + collagen_y_intercept)
                      << "\tEffective slope = " << kappa * eta_collagen
                      << "\tRest len = " << rest_len
                      << "\n";
        }*/ 
        return kappa * (eta_collagen*E + collagen_y_intercept);
    }
    else if (E > 0.0){
        /*if ((lag_mastr_idx % 2500) == 0){
            std::cout << "Exp.   (idx,nbr) = (" << lag_mastr_idx << ", " <<  lag_slave_idx << ")"
                      << "\tE = " << E
                      << "\tF = " << kappa * a * (exp(b*E) - 1)
                      << "\tEffective slope = " << kappa * a * b // taylor series coefficient on first term
                      << "\tRest len = " << rest_len
                      << "\n";
        }*/ 
        return kappa * a * (exp(b*E) - 1);
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





