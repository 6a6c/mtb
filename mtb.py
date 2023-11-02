import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.sparse as sp
from scipy.integrate import odeint, solve_ivp
from matplotlib import animation
from multiprocessing.pool import ThreadPool
import pandas as pd


def populateA(width, r, s, delta):
    n = width ** 2
    w = width
    
    P = sp.lil_matrix((n,n), dtype=np.float64)

    for i in range(n):
        
        #if(i % 1000 == 0): print(i)

        # we go:
        # p_00->00 p_10->00 p_20->00 ... p_01->00 p_11->00 p_21->00 ... p_02->00 p_12->00 p_22->00 ...
        # p_00->10 p_10->10 p_20->10 ... p_01->10 p_11->10 p_21->10 ...
        
        # plasmid count given by x % w
        # free count given by x / w (int)
        
        tot = 0
        # plasmid death to get to state i
        if(i + 1 < n): 
            P[i, i + 1] = ((i + 1) % w) * delta
            tot += (i % w) * r * (1 - s)
            
        # plasmid birth to get to state i
        if(i - 1 >= 0): 
            P[i, i - 1] = ((i - 1) % w) * r * (1 - s)
            tot += (i % w) * delta
            
        # free death to get to state i
        if(i + w < n): 
            P[i, i + w] = (int((i + w) / w)) * delta
            tot += int(i / w) * r + (i % w) * r * s
            
        # free birth to get to state i
        if(i - w >= 0): 
            P[i, i - w] = ((int((i - w) / w)) * r + (i % w) * r * s)
            tot += int(i / w) * delta
            
        # A_kk = - sum{l\neqk}(A_lk)
        P[i, i] = - tot
                
    #print(P)
    
    return P
    
    
    
def createp(width, starting_p, starting_f):
    p = np.zeros(width ** 2, dtype=np.float64)
    p[starting_f * width + starting_p] = 1
    
    return p


def plotp(p, fn):
    fig1 = plt.figure(dpi = 500)
    a = plt.imshow(scanline_column_vector(p))
    plt.colorbar(a)
    plt.title(fn)
    plt.savefig(fn)
        
def plott(tots, fn):
    global width
    fig1 = plt.figure(dpi = 500)
    
    x = range(width * 2)
    plt.plot(x, tots)
    plt.title(fn)
    plt.savefig(fn)

def plotr(inds, vals, fn):
    fig1 = plt.figure(dpi = 500)
    plt.scatter(inds, vals, s = 1)
    plt.axis([0, 1, 0, max(vals)])
    plt.title(fn)
    plt.savefig(fn) 

def scanline_column_vector(p):
    w = int(np.sqrt(np.size(p)))
    arr = []
    for i in range(w):
        l = p[(i*w):((i+1)*w)]
        arr.append(l)
    arr.reverse()
    
    arr = np.array(arr)
    return arr
    
def get_totals(p):
    w = int(np.sqrt(np.size(p)))
    arr = np.zeros(2 * w)
    for i in range(len(p)):
        arr[i % w + int(i / w)] += p[i]
    
    return arr

def get_ratios(p):
    w = int(np.sqrt(np.size(p)))
    vals = np.zeros(int(np.size(p)))
    inds = np.zeros(int(np.size(p)))
    for i in range(len(p)):
        if(i == 0): 
            vals[i] += p[i]
            inds[i] = 0
        else:
            vals[i] += p[i]
            inds[i] = (i % w) / (i / w + i % w)
            
    return vals, inds

def store_P(P, fn):
    df = pd.DataFrame(P)
    print(df)
    df.to_csv(fn)

def read_P(width, fn):
    P = np.zeros(width ** 2, dtype = np.float64)
    
    df = pd.read_csv(fn)
    
    for i in df.index:
        P[i] = df.iloc[i][1]
        
    return P

def read_P_expand(olw, new, fn):
    P = np.zeros(new ** 2, dtype = np.float64)
    
    df = pd.read_csv(fn)

    off = 0
    for i in df.index:
        if(i % olw == 0 and i != 0): off += new
        
        P[i + off] = df.iloc[i][1]


def get_avg_rat(P):
    w = int(np.sqrt(np.size(P)))
      
    ret = 0.0
    #tot = 0.0
    for i in range(1, len(P)):
        #if(i % w == 0): continue
        ret += P[i] * (i % w) / (int(i / w) + i % w)
        #tot += P[i]        

    return ret

def get_std_rat(P, avg):
    w = int(np.sqrt(np.size(P)))
    
    ret = 0.0
    #tot = 0.0
    for i in range(0, len(P)):
        if(i == 0):
            #continue
            ret += P[i] * (avg) ** 2
        else:
            ret += P[i] * (avg - ( (i % w) / (int(i / w) + i % w) ) ) ** 2
            #tot += P[i] 
    return ret

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
