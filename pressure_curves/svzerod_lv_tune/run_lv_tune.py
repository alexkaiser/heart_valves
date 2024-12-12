
import pysvzerod
from scipy.io import savemat
import json 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds

import pdb

MMHG_TO_CGS = 1333.22368


class optimizer_class():

    def __init__(self, 
                 fwd_sim_file, 
                 calibration_data_file, 
                 variables_to_opt, 
                 bounds_all,
                 cycle_duration=1.0,
                 maxit=1000,
                 plot_variables=None):

        with open(fwd_sim_file) as fname:
            self.fwd_sim_obj = json.load(fname)        

        with open(calibration_data_file) as fname:
            self.calibration_data = json.load(fname)                

        self.solver = pysvzerod.Solver(self.fwd_sim_obj)

        self.variables_to_opt = variables_to_opt

        if plot_variables is not None:
            self.plot_variables = plot_variables
        else:
            self.plot_variables = self.variables_to_opt

        self.solver_run = False

        self.cycle_duration = cycle_duration
        self.maxit = maxit

        # dictionary of all values used in optimization 
        self.values = {}
        for var_name in variables_to_opt:
            self.values[var_name] = None

        self.bounds_all = bounds_all



    def run_and_update(self):

        self.solver.run()
        for var_name in self.variables_to_opt:
            self.values[var_name] = np.array(self.solver.get_single_result(var_name))


    def objective_function(self):

        objective_value = 0.0

        for var_name in self.variables_to_opt:

            # time values for one cycle normalized for solver data 
            n_time_values = len(self.values[var_name])
            query_pts = np.linspace(0.0, self.cycle_duration, n_time_values)

            # time values for one cycle normalized for observation data 
            n_times_observation = len(self.calibration_data['y'][var_name])
            observation_pts = np.linspace(0.0, self.cycle_duration, n_times_observation)

            # linearly interpolated obs data
            data_interpolated = np.interp(query_pts, observation_pts, np.array(self.calibration_data['y'][var_name]))

            diffs = abs(self.values[var_name] - data_interpolated)
            # diffs = abs(self.values[var_name] - np.array(calibration_data_file['y'][var_name]))

            # no integration weights needed because cancel 
            diffs_two_norm_squared = (diffs ** 2).sum()
            data_two_norm_squared = (data_interpolated ** 2).sum() 

            rel_error = diffs_two_norm_squared / data_two_norm_squared

            objective_value += rel_error 

            # objective_value += (diffs ** 2).sum()

            # plot
            # fig, ax = plt.subplots()

            # ax.plot(query_pts, self.values[var_name], 'x', markeredgewidth=2)
            # ax.plot(query_pts, data_interpolated, linewidth=2.0)
            # plt.show()

        # print("objective_value = ", objective_value)

        return objective_value



    def plot(self, title_str=''):

        print("self.plot_variables = ", self.plot_variables)

        fig, axs = plt.subplots(len(self.plot_variables))
        fig.suptitle(title_str)

        for idx,var_name in enumerate(self.plot_variables):

            if "pressure" in var_name:
                scaling = 1.0/MMHG_TO_CGS
            else: 
                scaling = 1.0

            # time values for one cycle normalized for solver data 
            n_time_values = len(self.values[var_name])
            query_pts = np.linspace(0.0, self.cycle_duration, n_time_values)

            # time values for one cycle normalized for observation data 
            n_times_observation = len(self.calibration_data['y'][var_name])
            observation_pts = np.linspace(0.0, self.cycle_duration, n_times_observation)

            # linearly interpolated obs data
            data_interpolated = np.interp(query_pts, observation_pts, np.array(self.calibration_data['y'][var_name]))

            axs[idx].plot(query_pts, scaling * self.values[var_name], 'x', markeredgewidth=2)
            axs[idx].plot(query_pts, scaling * data_interpolated, linewidth=2.0)            
            axs[idx].set_title(var_name)

        plt.show(block=False)


    def output_mat(self, matfile_name):

        times = self.solver.get_times()
        q_aorta = self.solver.get_single_result("flow:vessel:OUTLET")
        p_aorta = self.solver.get_single_result("pressure:vessel:OUTLET")

        q_mitral = self.solver.get_single_result("flow:valve0:ventricle")
        p_lv_in = self.solver.get_single_result("pressure:valve0:ventricle")

        q_lv_out = self.solver.get_single_result("flow:ventricle:valve1")
        p_lv_out = self.solver.get_single_result("pressure:ventricle:valve1")

        q_aortic_valve = self.solver.get_single_result("flow:valve1:vessel")
        p_aortic_valve = self.solver.get_single_result("pressure:valve1:vessel")

        pressure_c_outlet = self.solver.get_single_result("pressure_c:OUTLET")

        volume_lv = self.solver.get_single_result("Vc:ventricle")

        # volume_lv_initial = volume_lv[0]
        # print("volume_lv_initial = ", volume_lv_initial)

        # for idx,t in enumerate(times):
        #     if t >= 0.2:
        #         print("at t = ", t, ", p_aortic_valve_mmHg = ", p_aorta[idx]/MMHG_TO_CGS)
        #         break 

        dic = {"times": times, 
               "q_aorta": q_aorta,
               "p_aorta": p_aorta,
               "q_mitral": q_mitral,
               "p_lv_in": p_lv_in,
               "q_lv_out": q_lv_out,
               "p_lv_out": p_lv_out,
               "q_aortic_valve": q_aortic_valve,
               "p_aortic_valve": p_aortic_valve,
               "pressure_c_outlet": pressure_c_outlet,
               "volume_lv": volume_lv
            }

        savemat("matfile_name", dic)


    def print_summary(self, targets):

        times = self.solver.get_times()
        p_aorta = self.solver.get_single_result("pressure:vessel:OUTLET")
        volume_lv = self.solver.get_single_result("Vc:ventricle")

        volume_lv_initial = volume_lv[0]
        print("volume_lv_initial = ", volume_lv_initial)

        for idx,t in enumerate(times):
            if t >= 0.2:
                print("at t = ", t, ", p_aortic_valve_mmHg = ", p_aorta[idx]/MMHG_TO_CGS)
                break 

        for target_name in targets:
            if target_name in self.fwd_sim_obj['chambers'][0]['values']:
                print(target_name, "= ", self.fwd_sim_obj['chambers'][0]['values'][target_name])
            elif target_name in self.fwd_sim_obj['boundary_conditions'][1]['bc_values']:
                print(target_name, "= ", self.fwd_sim_obj['boundary_conditions'][1]['bc_values'][target_name])



def objective_function_free_params(params, opt_instance, targets):

    # set new values & re-run 0d solver 
    for idx,target_name in enumerate(targets):
        if target_name in opt_instance.fwd_sim_obj['chambers'][0]['values']:
            opt_instance.fwd_sim_obj['chambers'][0]['values'][target_name] = params[idx]
        elif target_name in opt_instance.fwd_sim_obj['boundary_conditions'][1]['bc_values']:
            opt_instance.fwd_sim_obj['boundary_conditions'][1]['bc_values'][target_name] = params[idx]
        else: 
            raise ValueError("target must be in OUTLET boundary condition or chamber")

    opt_instance.solver = pysvzerod.Solver(opt_instance.fwd_sim_obj)

    try:
        opt_instance.run_and_update()
    except RuntimeError:
        return float('inf')

    return opt_instance.objective_function()    


def run_optimization(opt_instance, targets):

    x0 = []
    bounds = []    
    for idx,target_name in enumerate(targets):

        if target_name in opt_instance.fwd_sim_obj['chambers'][0]['values']:
            x0.append(opt_instance.fwd_sim_obj['chambers'][0]['values'][target_name])

        elif target_name in opt_instance.fwd_sim_obj['boundary_conditions'][1]['bc_values']:
            x0.append(opt_instance.fwd_sim_obj['boundary_conditions'][1]['bc_values'][target_name])
            # only work on outlet bc 
            assert opt_instance.fwd_sim_obj['boundary_conditions'][1]['bc_name'] == 'OUTLET'
        else: 
            raise ValueError("target must be in OUTLET boundary condition or chamber")

        bounds.append(opt_instance.bounds_all[target_name])
    
    # for target_name in targets:
    #     bounds.append(opt_instance.bounds_all[target_name])

    # x0 = [opt_instance.fwd_sim_obj['chambers'][0]['values']['Emax']]

    # print("opt_instance.maxit = ", opt_instance.maxit)

    result = minimize(objective_function_free_params, 
                      x0, 
                      args=(opt_instance, targets), 
                      bounds = bounds,
                      method='Nelder-Mead',
                      options={'maxiter': opt_instance.maxit}
                      )

    return result


def run_rcr_simple():

    solver = pysvzerod.Solver("pulsatileFlow_R_RCR.json")
    solver.run()
    results = solver.get_full_result()
    print(results)

    results.to_csv('pulsatileFlow_R_RCR_output.csv')


def run_valve_tanh():

    solver = pysvzerod.Solver("valve_tanh.json")
    solver.run()
    results = solver.get_full_result()
    print(results)
    results.to_csv('valve_tanh.csv')    

    times = solver.get_times()
    
    q_upstream = solver.get_single_result("flow:INLET:upstream_vessel")
    p_upstream = solver.get_single_result("pressure:INLET:upstream_vessel")
    q_downstream = solver.get_single_result("flow:downstream_vessel:OUTLET")
    p_downstream = solver.get_single_result("pressure:downstream_vessel:OUTLET")

    q_valve_upstream = solver.get_single_result("flow:upstream_vessel:valve")
    p_valve_upstream = solver.get_single_result("pressure:upstream_vessel:valve")
    q_valve_downstream = solver.get_single_result("flow:valve:downstream_vessel")
    p_valve_downstream = solver.get_single_result("pressure:valve:downstream_vessel")

    dic = {"times": times, 
           "q_upstream": q_upstream, 
           "p_upstream": p_upstream,
           "q_downstream": q_downstream,
           "p_downstream": p_downstream,
           "q_valve_upstream": q_valve_upstream,
           "p_valve_upstream": p_valve_upstream,
           "q_valve_downstream": q_valve_downstream,
           "p_valve_downstream": p_valve_downstream
           }

    savemat("valve_tanh.mat", dic)


def run_calibrator_test():

    with open("steadyFlow_calibration.json") as fname:
        config = json.load(fname)
    cali = pysvzerod.calibrate(config)


def run_chamber():
    # solver = pysvzerod.Solver("chamber_elastance_valve_rcr.json")
    solver = pysvzerod.Solver("chamber_elastance_two_hill_valve_rcr.json")

    solver.run()

    results = solver.get_full_result()
    print(results)
    results.to_csv('chamber_elastance_valve_rcr.csv')


    times = solver.get_times()
    q_aorta = solver.get_single_result("flow:vessel:OUTLET")
    p_aorta = solver.get_single_result("pressure:vessel:OUTLET")

    q_mitral = solver.get_single_result("flow:valve0:ventricle")
    p_lv_in = solver.get_single_result("pressure:valve0:ventricle")

    q_lv_out = solver.get_single_result("flow:ventricle:valve1")
    p_lv_out = solver.get_single_result("pressure:ventricle:valve1")

    q_aortic_valve = solver.get_single_result("flow:valve1:vessel")
    p_aortic_valve = solver.get_single_result("pressure:valve1:vessel")

    pressure_c_outlet = solver.get_single_result("pressure_c:OUTLET")

    volume_lv = solver.get_single_result("Vc:ventricle")

    dic = {"times": times, 
           "q_aorta": q_aorta,
           "p_aorta": p_aorta,
           "q_mitral": q_mitral,
           "p_lv_in": p_lv_in,
           "q_lv_out": q_lv_out,
           "p_lv_out": p_lv_out,
           "q_aortic_valve": q_aortic_valve,
           "p_aortic_valve": p_aortic_valve,
           "pressure_c_outlet": pressure_c_outlet,
           "volume_lv": volume_lv
        }

    savemat("chamber_elastance_valve_rcr.mat", dic)

    
if __name__== "__main__":

    # example runs 
    # run_rcr_simple()
    # run_valve_tanh()
    # run_calibrator_test()
    # run_chamber()

    run_opt = True
    if run_opt:

        cycle_duration = 0.8
        maxit = 10000

        fwd_sim_file = "chamber_elastance_two_hill_valve_rcr.json"

        calibration_data_file = "chamber_elastance_valve_rcr_calibrate.json"

        var_to_opt = ["flow:ventricle:valve1", "pressure:ventricle:valve1", "pressure:vessel:OUTLET"]    
        plot_vars = ["flow:ventricle:valve1", "pressure:ventricle:valve1", "pressure:vessel:OUTLET"]    

        targets_all = ['Emax', 'Emin', 't_shift', 'tau_1', 'tau_2', 'm1', 'm2', 'C', 'Rd', 'Rp']

        bounds_all = {'Emax'    : [0.0, 1e4], 
                      'Emin'    : [0.0, 1e3], 
                      't_shift' : [0.0, cycle_duration], 
                      'tau_1'   : [0.0, cycle_duration], 
                      'tau_2'   : [0.0, cycle_duration], 
                      'm1'      : [0.0, 40.0], 
                      'm2'      : [0.0, 40.0],
                      'C'       : [0.0, 0.01],
                      'Rd'      : [0.0, 4000.0],
                      'Rp'      : [0.0, 500.0]
                      }



        optimizer = optimizer_class(fwd_sim_file, 
                                    calibration_data_file, 
                                    var_to_opt, 
                                    bounds_all,
                                    cycle_duration, 
                                    maxit,
                                    plot_vars)


        optimizer.run_and_update()
        obj_val_1 = optimizer.objective_function()
        optimizer.plot('before')
        print("obj_val_1 = ", obj_val_1)


        # one variable at a time     
        # targets = ['Emin', 'Emax', 'Vrd', 'Vrs', 't_shift']
        # targets = ['t_shift']
        # # targets = ['tau_1']
        # # targets = ['Emin']

        # result = run_optimization(optimizer, targets)
        # print("result = ", result)

        # targets = ['Emin']
        # result = run_optimization(optimizer, targets)    
        # print("result = ", result)    

        # targets = ['Emax']
        # result = run_optimization(optimizer, targets)    
        # print("result = ", result)    

        # targets = ['tau_1', 'tau_2', 'm1', 'm2']
        # result = run_optimization(optimizer, targets)    
        # print("result = ", result)    

        targets = ['Emax', 'Emin', 't_shift', 'tau_1', 'tau_2', 'm1', 'm2']
        optimizer.variables_to_opt = ["flow:ventricle:valve1", "pressure:ventricle:valve1"]
        result = run_optimization(optimizer, targets)    
        print("result = ", result)    
        optimizer.print_summary(targets_all)
        optimizer.plot('after chamber')

        targets = ['C', 'Rd', 'Rp']
        optimizer.variables_to_opt = ["flow:ventricle:valve1", "pressure:vessel:OUTLET"]    
        result = run_optimization(optimizer, targets)    
        print("result = ", result)    
        optimizer.print_summary(targets_all)        
        optimizer.plot('after rcr')

        optimizer.variables_to_opt = ["flow:ventricle:valve1", "pressure:ventricle:valve1", "pressure:vessel:OUTLET"]    
        result = run_optimization(optimizer, targets_all)    
        print("result = ", result)    
        optimizer.print_summary(targets_all)                


        # print("resulting Emax = ", optimizer.fwd_sim_obj['chambers'][0]['values']['Emax'])
        # optimizer.run_and_update()
        # for target_name in targets:
        #     if target_name in optimizer.fwd_sim_obj['chambers'][0]['values']:
        #         print(target_name, "= ", optimizer.fwd_sim_obj['chambers'][0]['values'][target_name])
        #     elif target_name in optimizer.fwd_sim_obj['boundary_conditions'][1]['bc_values']:
        #         print(target_name, "= ", optimizer.fwd_sim_obj['boundary_conditions'][1]['bc_values'][target_name])

        optimizer.plot('after')

        optimizer.output_mat("chamber_elastance_two_hill_valve_rcr_optimized_results.mat")

        plt.show()

