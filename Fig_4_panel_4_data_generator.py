import xgi
import random
import math as math
import pandas as pd
import numpy as np
import scipy as sp
import time
import matplotlib.pyplot as plt
from CNO import CNO
from simulation import simulation
########## primary network ############
#####################################################################################
n = 1000                                              # network size
m = 2                                                # network order: pairwise
k = {i: 30 for i in range(n)}    # generates degree distribution 
H = xgi.uniform_hypergraph_configuration_model(k, m) # generates network
save_dict = {'network': False,
            'avg_opinion': True, 
            'frac': True
            }
gamma = 0.2                                          # recovery rate
betabar = 0.77                                  # infection rate
K = 0.1                                              # opinion shift rate
opinions = [0 for x in range(1,n+1)]                 # generates initial opinion distribution
#####################################################################################
################ simulation parameter dict ###############
#####################################################################################
simdict = {'network': H,
           'opinions': opinions,
           'infectionrate': betabar,
           'recoveryrate': gamma,
           'opinionshift': K,
           'ep': 2,
           'reinfectiontime': 0,
           'timescale': 0.25,
           'update_size_p1': math.ceil(0.3*n),
           'alt_update_size_p1': 60,
           'update_size_n1': math.ceil(0.3*n),
           'alt_update_size_n1': 60,
           'event_size':  math.ceil((4/15)*n),
           'save_dict': save_dict,
           'savefreq':  1,
           'exittime':  1500}
#####################################################################################

############### variable sweep ##############
#####################################################################################
print('Simulation for panel 4 has begun')
M = 3 # number of simulations per point
epvec = np.linspace(0.001,10,35)
xvec = np.linspace(-1,1,35)
samplevec = np.linspace(0,0.05,M)
sweepvar = np.zeros( (len(xvec), len(epvec)) )
sp1 = np.zeros(M*M)
sn1 = np.zeros(M*M)
i = 0
j = 0
for ep in epvec:
    i = 0
    for opinion in xvec:
        simdict['ep'] = ep
        simdict['opinions'] = [opinion for x in range(1,n+1)]
        q = 0
        for q1 in range(M):
            for q2 in range(M):
                simdict['update_size_p1'] = math.floor(n*samplevec[q1])
                simdict['update_size_n1'] = math.floor(n*samplevec[q2])
                sim = simulation(simdict)
                sim.main('single_infection')
                sp1[q] = sim.contagion_network.get_frac_state(1)
                sn1[q] = sim.contagion_network.get_frac_state(-1)
                q = q + 1
        if np.mean(sp1) > np.mean(sn1):
            sweepvar[i,j] = 1
        elif np.mean(sn1) > np.mean(sp1):
            sweepvar[i,j] = -1
        else:
            sweepvar[i,j] = 0
        i = i + 1
    j = j + 1

titlestring = 'sim_stability_sweep_P4.txt'
np.savetxt(titlestring, sweepvar)
print('Simulation for panel 4 has ended')
#####################################################################################