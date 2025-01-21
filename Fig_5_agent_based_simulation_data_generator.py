#################################################################################################################################################################################################################
# Copyright 2023, 2024, 2025 Corbit R. Sampson
#################################################################################################################################################################################################################
# This file is part of the Competing_Social_Contagions repository.
#
# Competing_Social_Contagions repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# Competing_Social_Contagions repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with the Competing_Social_Contagions repository. If not, see <https://www.gnu.org/licenses/>. 
#################################################################################################################################################################################################################


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
k = {i:7 for i in range(n)}    # generates degree distribution 
H = xgi.uniform_hypergraph_configuration_model(k, m) # generates network
save_dict = {'network': False,
            'avg_opinion': True, 
            'frac': True
            }
gamma = 0.2                                          # recovery rate
betabar = 0.584                                  # infection rate
K = 0.1001                                              # opinion shift rate
opinions = [0.3 for x in range(1,n+1)]                 # generates initial opinion distribution
#####################################################################################
################ simulation parameter dict ###############
#####################################################################################
simdict = {'network': H,
           'opinions': opinions,
           'infectionrate': betabar,
           'recoveryrate': gamma,
           'opinionshift': K,
           'ep': 1,
           'reinfectiontime': 0,
           'timescale': 0.25,
           'update_size_p1': math.ceil(0.5*n),
           'alt_update_size_p1': math.ceil(0.5*n),
           'update_size_n1': math.ceil(0.3*n),
           'alt_update_size_n1': math.ceil(0.3*n),
           'event_size':  math.ceil((4/15)*n),
           'save_dict': save_dict,
           'savefreq':  1,
           'exittime':  6000}
#####################################################################################

################ generate and run simulation ##################
#####################################################################################
sim = simulation(simdict,filename_ext = 'fig1data')
sim.main('single_infection')
#####################################################################################