# This file gets KS statistics for the combination/stochastic simulations
# against the ME solution.

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.integrate import odeint, solve_ivp
from multiprocessing.pool import ThreadPool
import pandas as pd
from numpy.random import choice
from scipy.stats import kstest, ks_2samp
import mtb

import time

def deriv(t, y):
    p = A * y
    return p


"""
Set these to be the r and delta that you want for comnination model
You can also specifiy a piecewise or other function.
"""
def delta(t):
    return .64
   
    
def r(t):
    return .24



def dNdt(t, N):
    return (r(t) - delta(t)) * N


def dfdt(t, f):
    s = .1
    return -s * r(t) * f

def getCdfN(arr):
    global cdfN
    cdf = []
    for i in arr:
        if(i > len(cdfN) ):
            cdf.append(1)
        else:
            cdf.append(cdfN[i])
    
    return cdf

if __name__ == "__main__":
    global cdfN
    directory = "data_r.84"
    width = 1000
    #A = mtb.populateA(width, r, s, delta)
    #P = mtb.createp(width, start_P, start_F)


    stochdf = pd.read_csv("gil_sims/simulations-P0-1-s-0.1-r-0.84-delta-0.44.csv")

    f = open("r84_ks_com.csv", "w")
    f.write("start,end,N_ks,N_p,f_ks,f_p\n")
    f.close()    
 
    f = open("r84_ks_stoch.csv", "w")
    f.write("day,N_ks,N_p,f_ks,f_p\n")
    f.close()    

    for start in range(1, 28):
        

        """ KS tests of current day stochastic simulations against current day ME soltuion """        
        startP = mtb.read_P(width, directory + f"/day_{start:.2f}_P.csv")
        
        state_Ns = [ ( (i % width) + (i // width) ) for i in range(len(startP)) ]
        state_fs = [0] + [ ( (i % width) / ((i % width) + (i // width)) ) for i in range(1, len(startP) ) ]
        

        me_N_choice = choice(state_Ns, 10000, p = abs(startP)/(sum(abs(startP))))
        me_f_choice = choice(state_fs, 10000, p = abs(startP)/(sum(abs(startP))))

        # This is why R is the worst. Who the hell makes csv's like this what the fuck.
        stoch_N_res = stochdf.loc[abs(stochdf["day"] - start) < .2]["total"].to_numpy()
        stoch_f_res = stochdf.loc[abs(stochdf["day"] - start) < .2 ]["fp"].to_numpy()


        stoch_N_res = np.concatenate( ( stoch_N_res, np.zeros( 1000 - len(stoch_N_res) ) ) )
        stoch_f_res = np.concatenate( ( stoch_f_res, np.zeros( 1000 - len(stoch_f_res) ) ) )
        
        print(len(stoch_N_res), stoch_N_res)

        
        Nks = ks_2samp(stoch_N_res, me_N_choice)
        print(Nks)
        fks = ks_2samp(stoch_f_res, me_f_choice)
        print(fks)

        with open("r84_ks_stoch.csv", "a") as f:
            s = str(start) + "," + str(Nks.statistic) + "," + str(Nks.pvalue) + "," + str(fks.statistic) + "," + str(fks.pvalue) + "\n"
            f.write(s)

    

        """ KS tests of 
            combination model with time switch at day start
            against ME solution at day end 
        """
        for end in range(start, 16):
            print(start, end)             
            endP = mtb.read_P(width, directory + f"/day_{end:.2f}_P.csv") 
            print(sum(endP), sum(abs(endP)))
             
            print("Ps read")
            print("states")    
            
            """
            Calculating the actual CDF for values has not been working. Ive tried literally everything
            It should work, have no idea wwhat im doing wrrong
            Instead, we're just gonna use two samples ks tests.
            These ones get the low values we expect from identical samples
            im so tired
         
            cdfN = {}; cdff = {}
            for i in range(len(endP)):
                if (i %width) + (i//width) in cdfN:
                    cdfN[ (i%width) + (i//width) ] += endP[i]
                else:
                    cdfN[(i%width)+(i//width)] = endP[i]
                if(i == 0):
                    cdff[0] = endP[i]        
                else:
                    if (i%width) / ((i % width) + (i // width)) in cdff:
                        cdff[(i%width) / ((i % width) + (i // width))] += endP[i] 
                    else:
                        cdff[(i%width) / ((i % width) + (i // width))] = endP[i] 

            cdfN = list(cdfN.values())
            cdff = list(cdff.values())
            print("N", len(cdfN), cdfN[-1], "f", len(cdff), cdff[-1])

            running = 0.0
            for i in range(len(cdfN)):
                running += cdfN[i]
                cdfN[i] = running
            running = 0.0
            for i in range(len(cdff)):
                running += cdff[i]
                cdff[i] = running
                 
            print("N", cdfN[-1], "f", cdff[-1])    
            """
            me_N_choice = choice(state_Ns, 10000, p = abs(endP)/(sum(abs(endP))))
            me_f_choice = choice(state_fs, 10000, p = abs(endP)/(sum(abs(endP))))

            #print("end choices")

            com_N_choice = choice(state_Ns, 10000, p = abs(startP)/(sum(abs(startP))))
            com_f_choice = choice(state_fs, 10000, p = abs(startP)/(sum(abs(startP))))
            print("start choices")

            com_N_res = np.zeros(10000)
            com_f_res = np.zeros(10000)
            for i in range(10000):
                if(start == end):
                    com_N_res[i] = com_N_choice[i]
                    com_f_res[i] = com_f_choice[i]
                    continue

                tmpN = solve_ivp(dNdt, [start, end], [ com_N_choice[i] ] , t_eval = [end] )
                    
                tmpf = solve_ivp(dfdt, [start, end], [ com_f_choice[i] ] , t_eval = [end])
                
                com_N_res[i] = tmpN.y
                
                com_f_res[i] = tmpf.y

            print("ivp solved")


            Nks = ks_2samp(com_N_res, me_N_choice)
            print(Nks)
            fks = ks_2samp(com_f_res, me_f_choice)
            print(fks)

            with open("r84_ks_com.csv", "a") as f:
                s = str(start) + "," + str(end) + "," + str(Nks.statistic) + "," + str(Nks.pvalue) + "," + str(fks.statistic) + "," + str(fks.pvalue) + "\n"
                f.write(s)

