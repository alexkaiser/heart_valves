import os, re
import numpy as np 
import matplotlib.pyplot as plt

import warnings

debug = True

def compute_rcr_parameters(R, P_min, P_max, ratio_prox_to_distal_resistors, decay_time, C_prefactor=1.0):

    tol = 1e-14

    # ratio of resistors is constant 
    # resistors sum to total resistance 
    R_distal = R / (1.0 + ratio_prox_to_distal_resistors)
    R_proximal = R - R_distal

    assert abs(R_distal + R_proximal - R) < tol 

    # timescale for pressure decrease during aortic valve closure 
    C = -C_prefactor * decay_time / (R_distal * np.log(P_min/P_max))

    return R_proximal, C, R_distal


if __name__== "__main__":

    MMHG_TO_CGS = 1333.22368

    # experimental mean pressure 
    p_patrunk_mean_mmHg = 38.617063756063736
    p_patrunk_mean = p_patrunk_mean_mmHg * MMHG_TO_CGS

    # downstream pressure is assumed to be ground, zero 

    Q_mean_total_l_per_min = 3.5 

    # total flow mean 
    Q_mean = Q_mean_total_l_per_min * 1000.0 / 60 

    # total downstream resistance 
    R_total_mmHg = p_patrunk_mean_mmHg / Q_mean
    R_total = p_patrunk_mean / Q_mean

    print("Q_mean = ", Q_mean)
    print("R_total_mmHg = ", R_total_mmHg)

    # steady flow sim numbers 
    # NOT CONVERGED! TEMP VALUES!! 
    warnings.warn("temporary values of steady sim being used")

    deltap_steady_mmHg = 1.0
    deltap_steady = deltap_steady_mmHg * MMHG_TO_CGS
    q_lpa_steady_flow_final = 35.475201349000002
    q_rpa_steady_flow_final = 23.474931795000000
    q_total_steady_flow = 58.944214744000000

    ratio_R_rpa_over_R_lpa = q_lpa_steady_flow_final / q_rpa_steady_flow_final

    print("ratio_R_rpa_over_R_lpa = ", ratio_R_rpa_over_R_lpa)

    # resistance of the lpa segment via steady sim 
    R_lpa = (1.0/ratio_R_rpa_over_R_lpa + 1.0) * deltap_steady / q_total_steady_flow
    R_lpa_mmHg = (1.0/ratio_R_rpa_over_R_lpa + 1.0) * deltap_steady_mmHg / q_total_steady_flow

    # right comes from the ratio 
    R_rpa = R_lpa * ratio_R_rpa_over_R_lpa
    R_rpa_mmHg = R_lpa_mmHg * ratio_R_rpa_over_R_lpa

    print("R_lpa = ", R_lpa, "R_lpa_mmHg = ", R_lpa_mmHg)
    print("R_rpa = ", R_rpa, "R_rpa_mmHg = ", R_rpa_mmHg)

    # from the total resistance formula 
    R_outlet_lpa      = 2.0 * R_total - R_lpa
    R_outlet_lpa_mmHg = 2.0 * R_total_mmHg - R_lpa_mmHg

    # from requirement that total resitance on each side is equal 
    R_outlet_rpa = R_lpa - R_rpa + R_outlet_lpa
    R_outlet_rpa_mmHg = R_lpa_mmHg - R_rpa_mmHg + R_outlet_lpa_mmHg

    print("R_outlet_lpa = ", R_outlet_lpa, "R_outlet_lpa_mmHg = ", R_outlet_lpa_mmHg)
    print("R_outlet_rpa = ", R_outlet_rpa, "R_outlet_rpa_mmHg = ", R_outlet_rpa_mmHg)

    decay_time = .4 

    # values from the experimental trace 
    P_max = 42 * MMHG_TO_CGS # interpolating the decay by eye to the middle of the oscillation 
                             # in the experimental trace of PA pressure 
    P_min = 31.25 * MMHG_TO_CGS 

    ratio_prox_to_distal_resistors = 0.1

    C_prefactor = 1.0

    R_p_rpa, C_rpa, R_d_rpa = compute_rcr_parameters(R_outlet_rpa, P_min, P_max, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

    R_p_lpa, C_lpa, R_d_lpa = compute_rcr_parameters(R_outlet_lpa, P_min, P_max, ratio_prox_to_distal_resistors, decay_time, C_prefactor)

    print("right_pa_R_proximal = %.14f" % R_p_rpa)
    print("right_pa_C = %.14f" % C_rpa)
    print("right_pa_R_distal = %.14f" % R_d_rpa)
    print("left_pa_R_proximal = %.14f" % R_p_lpa)
    print("left_pa_C = %.14f" % C_lpa)
    print("left_pa_R_distal = %.14f" % R_d_lpa)

