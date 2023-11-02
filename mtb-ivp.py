import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse as sp
from scipy.integrate import odeint, solve_ivp
from matplotlib import animation
from multiprocessing.pool import ThreadPool
import pandas as pd

import mtb


def deriv(t, y):
    p = A * y
    return p

if __name__ == "__main__":

    width = 1000
    r = 1.24
    delta = .84
    s = .1
    dt = .001
    start_P = 1
    start_F = 0

    A = mtb.populateA(width, r, s, delta)
    P = mtb.createp(width, start_P, start_F)


    for i in range(0, 112):
        tmp = solve_ivp(deriv, [i * .25, (i + 1) * .25], P, t_eval = [(i + 1) * .25])
        P = tmp.y[:,0]
        
        print(i)
        print(np.sum(P))
       
        mtb.store_P(P, "data_r1.24/day_{:.2f}_P.csv".format((i + 1) * .25))
