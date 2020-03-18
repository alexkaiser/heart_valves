import os, re
import numpy as np 
import matplotlib.pyplot as plt

debug = True

def compute_rcr_parameters(P_mean, Q_min, Q_max, Q_mean, dp_dt, frac_proximal, frac_distal):

    R_total = P_mean / Q_mean 

    C = (Q_max - Q_min) / dp_dt

    tol = 1e-14
    assert abs(frac_proximal + frac_distal - 1) < tol 

    R_p = frac_proximal * R_total
    R_d = frac_distal * R_total

    return R_p, C, R_d, R_total


if __name__== "__main__":

    MMHG_TO_CGS = 1333.22368

    P_mean = 1.0433437500000039e+02 * .9
    dp_dt  = 1300.0 * .25
    Q_min  = 0
    Q_max  = 500 
    # 5.6 L/min
    Q_mean = 5600.0 / 60 # ml/s 
    frac_proximal = 0.05
    frac_distal = 0.95

    R_p_mmHg, C_mmHg, R_d_mmHg, R_total_mmHg = compute_rcr_parameters(P_mean, Q_min, Q_max, Q_mean, dp_dt, frac_proximal, frac_distal)

    name = 'aorta'

    print "Values mmHg"
    print name, ",\t", R_p_mmHg, ",\t", C_mmHg, ",\t", R_d_mmHg, ",\t", R_total_mmHg
    print "\n\n\n"

    P_mean *= MMHG_TO_CGS
    dp_dt *= MMHG_TO_CGS
    R_p, C, R_d, R_total = compute_rcr_parameters(P_mean, Q_min, Q_max, Q_mean, dp_dt, frac_proximal, frac_distal)

    print "Values CGS"
    print name, ",\t", R_p, ",\t", C, ",\t", R_d, ",\t", R_total

    print "R_proximal = ", R_p
    print "C = ", C
    print "R_distal = ", R_d