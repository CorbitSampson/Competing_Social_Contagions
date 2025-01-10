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
#You should have received a copy of the GNU General Public License along with Foobar. If not, see <https://www.gnu.org/licenses/>. 
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
from opinion_generator import *
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
betabar = 0.6                                # infection rate
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
           'exittime':  3000}
#####################################################################################

############### variable sweep ##############
#####################################################################################
print('Simulation for panel 4 has begun')
M = 2 # number of simulations per point
epvec = np.linspace(0.001,10,35)
xvec = np.linspace(-1,1,35)
samplevec = np.linspace(0,0.05,M)
sweepvar = np.zeros( (len(xvec), len(epvec)) )
sp1 = np.zeros(M*M)
sn1 = np.zeros(M*M)
change_op = np.zeros(M*M)
change_sn1 = np.zeros(M*M)
change_sp1 = np.zeros(M*M)
change0max = np.zeros( (len(xvec), len(epvec)) )
changepmax = np.zeros( (len(xvec), len(epvec)) )
changenmax = np.zeros( (len(xvec), len(epvec)) )
change0avg = np.zeros( (len(xvec), len(epvec)) )
changepavg = np.zeros( (len(xvec), len(epvec)) )
changenavg = np.zeros( (len(xvec), len(epvec)) )

i = 0
j = 0
for ep in epvec:
    i = 0
    for opinion in xvec:
        simdict['ep'] = ep
        simdict['opinions'] = normal_wrapped(opinion, 0.5, n)
        q = 0
        for q1 in range(M):
            for q2 in range(M):
                simdict['update_size_p1'] = math.floor(n*samplevec[q1])
                simdict['update_size_n1'] = math.floor(n*samplevec[q2])
                sim = simulation(simdict, filename_ext = 'P4')
                sim.main('single_infection')
                sp1[q] = sim.contagion_network.get_frac_state(1)
                sn1[q] = sim.contagion_network.get_frac_state(-1)
                change_op[q] = np.abs(sim.avg_opinion[-1] - np.mean(sim.avg_opinion[-200:-1]))
                change_sn1[q] = np.abs(sim.frac_sn1[-1] - np.mean(sim.frac_sn1[-200:-1]))
                change_sp1[q] = np.abs(sim.frac_sp1[-1] - np.mean(sim.frac_sp1[-200:-1]))
                q = q + 1
        if np.mean(sp1) > np.mean(sn1):
            sweepvar[i,j] = 1
        elif np.mean(sn1) > np.mean(sp1):
            sweepvar[i,j] = -1
        else:
            sweepvar[i,j] = 0
        
        change0max[i,j] = np.max(change_op)
        changepmax[i,j] = np.max(change_sp1)
        changenmax[i,j] = np.max(change_sn1)
        change0avg[i,j] = np.mean(change_op)
        changepavg[i,j] = np.mean(change_sp1)
        changenavg[i,j] = np.mean(change_sn1)
        
        i = i + 1
    j = j + 1

titlestring = 'het_sim_stability_sweep_P4.txt'
np.savetxt(titlestring, sweepvar)
print('Simulation for panel 4 has ended')
print('Max change after 3000 time steps is:')
#np.savetxt('max_change_0_P4_H.txt',change0max)
#np.savetxt('max_change_p_P4_H.txt',changepmax)
#np.savetxt('max_change_n_P4_H.txt',changenmax)
#np.savetxt('avg_change_0_P4_H.txt',change0avg)
#np.savetxt('avg_change_0_P4_H.txt',changepavg)
#np.savetxt('avg_change_0_P4_H.txt',changenavg)
print(np.max(np.max(change0max)))
print(np.max(np.max(changepmax)))
print(np.max(np.max(changenmax)))
print(np.mean(np.mean(change0avg)))
print(np.mean(np.mean(changepavg)))
print(np.mean(np.mean(changenavg)))
#####################################################################################