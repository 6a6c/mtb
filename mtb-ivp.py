# This file solves the IVP for the master equation on the CPU.
# It creates the transition matrix A and proability vector P and then
# uses solve_ivp to solve for every .25 days

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse as sp
from scipy.integrate import odeint, solve_ivp
from matplotlib import animation
from multiprocessing.pool import ThreadPool
import pandas as pd

import mtb

import time

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

    A = mtb.populateA(width, r, s, delta)
    P = mtb.createp(width, start_P, start_F)
    #P = mtb.read_P(width, "data/day_28.00_P.csv")

    start = time.time()

    for i in range(0, 112):
        tmp = solve_ivp(deriv, [i * .25, (i + 1) * .25], P, t_eval = [(i + 1) * .25])
        P = tmp.y[:,0]
        
        print(i)
        print(np.sum(P))
      
        # Uncomment and change if you'd like to store intermediary data 
        #mtb.store_P(P, "data_det/day_{:.2f}_P.csv".format((i + 1) * .25))

    end = time.time()
   
    print("r value",  r) 
    print("FINAL TIME IN SECONDS: ", end - start)
