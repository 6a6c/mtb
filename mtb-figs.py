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
    std = mtb.get_std_rat_noex(P, avg):
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
    


 
def plotall(P, i):

    vtmp, itmp = mtb.get_ratios(P)
    ttmp = mtb.get_totals(P)
    
    mtb.plotp(P, dirt + "/figs/day_{:.2f}_contour.png".format(i* .25), i *.25)
    mtb.plott(ttmp, dirt + "/figs/day_{:.2f}_totals.png".format(i*.25), i*.25)
    mtb.plotr(itmp, vtmp, dirt + "/figs/day_{:.2f}_ratios.png".format(i*.25), i *.25)  


##### BEGIN MAIN
if __name__ == "__main__":

    width = 1000
    r = 1.04
    delta = .64
    s = .1
    dt = .001
    start_P = 1
    start_F = 0
    dirt = "res.64"




    x = np.zeros(113)

    stds = np.zeros(113)
    avgs = np.zeros(113)
    avgs_noex = np.zeros(113)
    stds_noex = np.zeros(113)

    tas = np.zeros(113)
    tss = np.zeros(113)
    tasn = np.zeros(113)
    tssn = np.zeros(113)

    P = mtb.createp(width, start_P, start_F)

    getall(P, 0)
    plotall(P, 0)

    for i in range(1,113):
        P = mtb.read_P(width, "data_r.64/day_{:.2f}_P.csv".format(i * .25))
        getall(P, i)
        plotall(P, i)


    print(stds)
    print(avgs)
    print(stds_noex)
    print(avgs_noex)
    print(tas)
    print(tss)
    print(tasn)
    print(tssn)

    df = pd.DataFrame({'Day': x, 'AvgRat': avgs, 'StdRat': stds, 'AvgRatNoEx': avgs_noex, 'StdRatNoEx': stds_noex, 'AvgTot': tas, 'StdTot': tss, 'AvgTotNoEx': tasn, 'StdTotNoEx': tssn})
    df.to_csv(dirt + "/info.csv")


    fig1 = plt.figure(dpi = 500)
    plt.title("Avg Fraction of Plasmid-Bearing Cells per Day")
    plt.ylabel("Avg Fraction")
    plt.xlabel("Day")
    plt.plot(x, avgs, color = "red")
    plt.fill_between(x, avgs - stds, avgs + stds, color = "red", alpha = .4)
    plt.savefig(dirt + "/rat.png")
    plt.close()

    fig1 = plt.figure(dpi = 500)
    plt.title("Avg Fraction of Plasmid-Bearing Cells per Day - No Extinction")
    plt.ylabel("Avg Fraction")
    plt.xlabel("Day")
    plt.plot(x, avgs_noex, color = "red")
    plt.fill_between(x, avgs_noex - stds_noex, avgs_noex + stds_noex, color = "red", alpha = .4)
    plt.savefig(dirt + "/rat_noex.png")
    plt.close()

    fig1 = plt.figure(dpi = 500)
    plt.title("Avg Total Cell Count per Day")
    plt.ylabel("Avg Total")
    plt.xlabel("Day")
    plt.plot(x, tas, color = "red")
    plt.fill_between(x, tas - tss, tas + tss, color = "red", alpha = .4)
    plt.savefig(dirt + "/tot.png")
    plt.close()

    fig1 = plt.figure(dpi = 500)
    plt.title("Avg Total Cell Count per Day - No Extinction")
    plt.ylabel("Avg Total")
    plt.xlabel("Day")
    plt.plot(x, tasn, color = "red")
    plt.fill_between(x, tasn - tssn, tasn + tssn, color = "red", alpha = .4)
    plt.savefig(dirt + "/tot_noex.png")
    plt.close()

