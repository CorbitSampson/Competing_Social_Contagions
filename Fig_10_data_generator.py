import xgi
import random
import pandas as pd
import numpy as np
import scipy as sp
import time
import matplotlib.pyplot as plt
from reduced_meanfield_object import reduced_meanfield_object

out_param = [0.1,0.5]
other_param_vec = [0.3,0.9]
p = 0

initial_condition_infectionP = np.linspace(0,1,50)
initial_condition_infectionN = np.linspace(0,1,50)
T = np.zeros((len(initial_condition_infectionN), len(initial_condition_infectionP)))

tend = 3000

for initialP in out_param:
    m = (4/15)
    C = 0.1
    betabar = 0.6
    mk = 6
    gamma = 0.2
    ep = 0.01
    x = 0.1
    r0 = (m*(mk-1)*betabar)/(mk*gamma)
    K = 0.25
    timeconstant = 0.1
    
    for param in other_param_vec:
        j = 0
        for x in initial_condition_infectionP:
            i = 0
            for forcing in initial_condition_infectionN:
                r0 = param
                v0 = [initialP,0,x]
                system = reduced_meanfield_object(r0,K,ep,timeconstant, v0, 'constant_function')
                t,y = system.euler_solver([0,tend], 0.25,value=forcing)
                T[i,j] = y[2,-1]
                i = i + 1
            j = j + 1

        titlestring = f'numeric_constant_force_data_{p}.txt'
        np.savetxt(titlestring,T)
        p = p + 1