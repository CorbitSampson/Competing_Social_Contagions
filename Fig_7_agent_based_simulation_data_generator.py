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
n = 3000                                              # network size
m = 2                                                # network order: pairwise
k = {i:10 for i in range(n)}    # generates degree distribution 
H = xgi.uniform_hypergraph_configuration_model(k, m) # generates network
save_dict = {'network': False,
            'avg_opinion': True, 
            'frac': True
            }
gamma = 0.2                                          # recovery rate
betabar = 0.6                                  # infection rate
K = 0.1                                              # opinion shift rate
opinions = [0.3 for x in range(1,n+1)]                 # generates initial opinion distribution
#####################################################################################
################ simulation parameter dict ###############
#####################################################################################
simdict = {'network': H,
           'opinions': opinions,
           'infectionrate': betabar,
           'recoveryrate': gamma,
           'opinionshift': K,
           'ep': 10,
           'reinfectiontime': 0,
           'timescale': 0.05,
           'update_size_p1': math.ceil(0.2*n),
           'alt_update_size_p1': 0,
           'update_size_n1': math.ceil(0.5*n),
           'alt_update_size_n1': 0,
           'event_size':  math.ceil((8/15)*n),
           'save_dict': save_dict,
           'savefreq':  1,
           'exittime':  9000}
#####################################################################################

################ generate and run simulation ##################
#####################################################################################
sim = simulation(simdict,filename_ext = 'fig2data')
sim.main('single_infection')
#####################################################################################