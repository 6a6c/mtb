# This file is set up to obtain various statistics like average, variance, 
# and percentiles from the ME solutions. Currently, its set up to create
# csv files for the percentiles of cell counts/plasmid fraction.
# It's pretty slow, mainly due to the cost of reading the solution files.
# Definitely comment out any calls/stats you don't wanna see

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse as sp
from scipy.integrate import odeint, solve_ivp
from matplotlib import animation
from multiprocessing.pool import ThreadPool
import pandas as pd

import mtb

def getall(P, i):
    x[i] = (i * .25)
    avg = mtb.get_avg_rat(P)
    std = mtb.get_std_rat(P, avg)
    avgs[i] = avg
    stds[i] = std
    print(i, "rat_ex", avg, std, end = " ")   
   
    avg = mtb.get_avg_rat_noex(P)
    std = mtb.get_std_rat_noex(P, avg)
    avgs_noex[i] = avg
    stds_noex[i] = std
    print("rat_noex", avg, std, end = " ")

    avg = mtb.get_avg_tot(P)
    std = mtb.get_std_tot(P, avg)
    tas[i] = avg
    tss[i] = std
    print("tot_ex", avg, std, end = " ")


    avg = mtb.get_avg_tot_noex(P)
    std = mtb.get_std_tot_noex(P, avg)
    tasn[i] = avg
    tssn[i] = std
    print("tot_noex", avg, std, end = "\n")
    


def get_log_avg_rat_noex(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    tot = 0.0
    for i in range(1, len(P)):
        #if(i % w == 0): continue
        if(100 * (i % w) / (int(i/w) + i %w) <= .1):
            continue
            #ret += P[i] * np.log10(.1)
        else:
            ret += P[i] * np.log10( 100 * (i % w) / (int(i / w) + i % w))
            tot += P[i]        

    return ret / tot



def get_log_avg_rat(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    #tot = 0.0
    for i in range(1, len(P)):
        #if(i % w == 0): continue
        if(100 * (i % w) / (int(i/w) + i %w) <= .1):
            ret += P[i] * np.log10(.1)
        else:
            ret += P[i] * np.log10( 100 * (i % w) / (int(i / w) + i % w))
        #tot += P[i]        

    return ret

def get_log_avg_tot(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    tot = 0.0
    for i in range(1, len(P)):
        ret += P[i] * np.log10(((i/w) + (i%w)))
        tot += P[i]        

    return ret

def get_log_avg_tot_noex(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    tot = 0.0
    for i in range(1, len(P)):
            
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            continue

        ret += P[i] * np.log10((i/w) + (i%w))
        tot += P[i]        

    return ret / tot


def get_tiles_rat_noex(P):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
    ma = 1

    for i in percorder:
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            ma -= P[i]


    for i in percorder:
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            continue
        tot += P[i]
        if((tot/ma) >= .05 and tile5 == -2):
            if(i == 0 or (i%w)/(int(i/w)+i%w) < .001 ):
                print("error somehow")
                tile5 = -1
            else:
                tile5 = np.log10( 100 * (i % w) / (int(i/w) + i%w))
        if((tot/ma) >= .95 and tile95 == -2):
            tile95 = np.log10( 100 * (i % w) / (int(i/w) + i%w))
            break

    return tile5, tile95


def get_tiles_rat(P):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
    for i in percorder:
        tot += P[i]
        if(tot >= .05 and tile5 == -2):
            if(i == 0 or (i%w)/(int(i/w)+i%w) < .001 ):
                tile5 = -1
            else:
                tile5 = np.log10( 100 * (i % w) / (int(i/w) + i%w))
        if(tot >= .95 and tile95 == -2):
            tile95 = np.log10( 100 * (i % w) / (int(i/w) + i%w))
            break

    return tile5, tile95

def get_tiles_tot_noex(P):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
    
    ma = 1
    for i in percorder:
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            ma -= P[i]

    for i in totorder:
        
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            continue
        tot += P[i]
        if((tot/ma) >= .05 and tile5 == -2):
                tile5 = np.log10((i%w) +(i/w))
        if((tot/ma) >= .95 and tile95 == -2):
            tile95 = np.log10((i%w) + (i/w))
            break

    return tile5, tile95


def get_tiles_tot(P):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
    for i in totorder:
        tot += P[i]
        if(tot >= .05 and tile5 == -2):    
            tile5 = np.log10((i/w) + (i%w))
        if(tot >= .95 and tile95 == -2):
            tile95 = np.log10((i/w) + (i%w))
            break

    return tile5, tile95

def get_avg_rat_noex(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    tot = 0.0
    for i in range(1, len(P)):
        if(i == 0): continue
        ret += P[i] * (i % w) / (int(i / w) + i % w)
        tot += P[i]        

    return ret / tot

def get_std_rat_noex(P, avg):
    w = int(np.sqrt(np.size(P)))
    
    ret = 0.0
    tot = 0.0
    for i in range(0, len(P)):
        if(i == 0):
            continue
            #ret += P[i] * (avg) ** 2
        else:
            ret += P[i] * (avg - ( (i % w) / (int(i / w) + i % w) ) ) ** 2

def get_tiles_tot_noex(P, percorder, totorder, desired_tots):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
     
    ma = 1
    """
    for i in percorder:
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            ma -= P[i]
    """

    
    #  Determines the total cell count percentiles
    ret = []
    working = desired_tots.pop(0)
    for i in totorder:
    
        if(i == 0):
            continue       
 
        tot += P[i]
        if((tot/ma) >= working):
                # Please note that we are going to take the first value
                # higher than the desired percentile as the percentile
                # It could be much higher, so this might cause problems
                ret.append(np.log10((i%w) +(i/w)))
                if(len(desired_tots) == 0):
                    break
                else:
                    working = desired_tots.pop(0)

    return ret

def get_tiles_rat_noex(P, percorder, totorder, desired_tots):
    w = int(np.sqrt(np.size(P)))
   
    tile5 = -2
    tile95 = -2
    tot = 0.0
    
    
    # Prunes out things beneath Limit of Detection
    ma = 1
    for i in percorder:
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            ma -= P[i]

    #  Determines the total cell count percentiles
    ret = []
    working = desired_tots.pop(0)
    for i in percorder:
        
        if(i == 0 or (i%w)/(int(i/w)+i%w) <= .001 ):
            continue
        tot += P[i]
        if((tot/ma) >= working):
                # Please note that we are going to take the first value
                # higher than the desired percentile as the percentile
                # It could be much higher, so this might cause problems
                ret.append(100 *  (i%w) / ((i%w) +(i/w)))
                if(len(desired_tots) == 0):
                    break
                else:
                    working = desired_tots.pop(0)

    return ret





##### BEGIN MAIN
if __name__ == "__main__":

    width = 1000
    r = .64         # modify this to be the r value you want to calc
    delta = .64
    s = .1
    dt = .001
    start_P = 1
    start_F = 0
    dirt = "res.64" # point this to the directory you have the ME
                    # solution files in.

    global percorder, totorder
    X1 = range(0, width ** 2)
    X2 = []
    X3 = []
    for i in X1:
        if(i == 0): 
            X2.append(0)
            X3.append(0)
        else:
            X2.append((i % width) / (int(i/width) + i %width))
            X3.append((i/width) + (i%width))

    percorder = [a for _,a in sorted(zip(X2, X1)) ] 
    totorder = [a for _,a in sorted(zip(X3, X1)) ]
    
    desired_tots = [ i/100 for i in range(1, 100) ]
    rs = [.64]
    
    dfr = pd.DataFrame(columns = desired_tots) 
    dft = pd.DataFrame(columns = desired_tots)

    for day in range(1,28):
        
        P = mtb.read_P(width, "data_r" + f"{r:.2f}".lstrip('0') + f"/day_{day}.00_P.csv")
        
        desired_tots = [ i/100 for i in range(1, 100) ]
        t = get_tiles_tot_noex(P, percorder, totorder, desired_tots)
        
        desired_tots = [ i/100 for i in range(1, 100) ]
        rat = get_tiles_rat_noex(P, percorder, totorder, desired_tots)
        
        print(t)
        print(rat)
        
        dfr.loc[day] = rat
        dft.loc[day] = t

        dfr.to_csv("rat64tiles.csv")
        dft.to_csv("tot64tiles.csv")

        
        ar = get_log_avg_tot_noex(P)
        print(f"Day {day} {r} log10 avg tot -NOEX: ", ar)

        ar = get_log_avg_rat_noex(P)
        print(f"Day {day} {r} log10 avg ratio - NOEX: ", ar)




