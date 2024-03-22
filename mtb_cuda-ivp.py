# This file is the same as mtb-ivp.py, but uses cupy for hardware acceleration
# Please clone this repo (in addition to setting up cupy) for the ivp solver:
# https://github.com/yor-dev/cupy_ivp/tree/master

import cupy as cp
import cupyx.scipy as csp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse as sp
from scipy.integrate import solve_ivp
from matplotlib import animation
from multiprocessing.pool import ThreadPool
import pandas as pd

import mtb_cuda
from cupy_ivp import solve_ivp, RK45

import time
import parser

def deriv(t, y):
    p = A * y
    return p

if __name__ == "__main__":

    width = 1000
    r = 1.04
    delta = .64
    s = .1
    dt = .001
    start_P = 1
    start_F = 0

    A = mtb_cuda.populateA(width, r, s, delta)
    P = mtb_cuda.createp(width, start_P, start_F)
    #P = mtb.read_P(width, "data/day_28.00_P.csv")

    start = time.time()

    for i in range(0, 112):
        tmp = solve_ivp(deriv, [i * .25, (i + 1) * .25], P, t_eval = [(i + 1) * .25])
        y = tmp.y[:,0]
        P = y       
 
        print(i)
        print(cp.sum(P))
    
        # Uncomment and change if you'd like to store intermediary data    
        #if(i % 4 == 3):
        #    mtb_cuda.store_P(P, "data_2000_r1.04/day_{:.2f}_P.csv".format((i + 1) * .25))

    
    end = time.time()

    print("r value: ", r)
    print("FINAL TIME IN SECONDS: ", end - start)


