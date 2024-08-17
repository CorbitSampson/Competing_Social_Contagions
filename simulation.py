# this file needs documentation


import xgi
import random
import pandas as pd
import numpy as np
import scipy as sp
import time
import json as json
import matplotlib.pyplot as plt
from CNO import CNO

class simulation:
    # argument dict = network, opinions, infectionrate, recoveryrate, timescale = 0.25, update_size_p1, update_szie_n1, savefreq, exittime = 1000
    def __init__(self, argument_dict,filename_ext = None):
        
        self.filename_ext = filename_ext
        self.update_size_p1 = argument_dict['update_size_p1']
        self.alt_update_size_p1 = argument_dict['alt_update_size_p1']
        self.update_size_n1 = argument_dict['update_size_n1']
        self.alt_update_size_n1 = argument_dict['alt_update_size_n1']
        self.spread_event_size = argument_dict['event_size']
        self.save_dict = argument_dict['save_dict']
        self.savefreq = argument_dict['savefreq']
        self.exittime = argument_dict['exittime']
        self.reinfectiontime = argument_dict['reinfectiontime']
        self.contagion_network = CNO(argument_dict['network'], 
                                argument_dict['opinions'],
                                argument_dict['infectionrate'],
                                argument_dict['recoveryrate'],
                                argument_dict['opinionshift'],
                                argument_dict['ep'],
                                argument_dict['timescale'])
        
        self.avg_opinion = []
        self.frac_sp1 = []
        self.frac_sn1 = []
        self.frac_s0  = []
        self.simtime = 0
        
    def time_step(self):
        spreaders = self.contagion_network.get_spreaders(self.spread_event_size)
        self.contagion_network.spread_event(spreaders)
        self.contagion_network.heal_event(self.contagion_network.network.nodes.filterby_attr('state', 0, mode='neq'))
        self.simtime = self.simtime + 1
        
    def save_network(self):
        if self.save_dict['network']:
            if self.simtime % self.savefreq == 0:
                titlestring = f'network_attime_{self.simtime}.json'
                xgi.write_json(self.contagion_network.network, titlestring)

    def exit_sim(self):
        if self.save_dict['avg_opinion']:
            if self.filename_ext == None:
                titlestring = 'avg_opinion.json'
            else:
                titlestring = 'avg_opinion_' + self.filename_ext + '.json'
            self.save_list(titlestring, self.avg_opinion)

        if self.save_dict['frac']:
            if self.filename_ext == None:
                titlestring1 = 'frac_sp1.json'
                titlestring2 = 'frac_sn1.json'
                titlestring3 = 'frac_s0.json'
            else:
                titlestring1 = 'frac_sp1_' + self.filename_ext + '.json' 
                titlestring2 = 'frac_sn1_' + self.filename_ext + '.json'
                titlestring3 = 'frac_s0_'  + self.filename_ext + '.json' 

            self.save_list(titlestring1, self.frac_sp1)
            self.save_list(titlestring2, self.frac_sn1)
            self.save_list(titlestring3, self.frac_s0)            

    def save_list(self, titlestring, list):
        with open(titlestring, 'w') as fp:
            json.dump(list, fp)

    def start_infection(self, first_infection = True):
        if first_infection:
            self.contagion_network.seed_contagion(self.update_size_p1, self.update_size_n1)
        else:
            self.contagion_network.seed_contagion(self.alt_update_size_p1, self.alt_update_size_n1, first_infection=first_infection)
    
    def store_states(self):
        if self.save_dict['avg_opinion']:
            self.store_avg_opinion()

        if self.save_dict['frac']:
            self.store_frac()

    def store_frac(self):
        self.frac_sp1.append(self.contagion_network.get_frac_state(1))
        self.frac_sn1.append(self.contagion_network.get_frac_state(-1))
        self.frac_s0.append(self.contagion_network.get_frac_state(0))

    def store_avg_opinion(self):
        self.avg_opinion.append(self.contagion_network.get_opinion_avg())

    def main(self, infection_mode):
        
        if infection_mode == 'reset_infection':
            self.start_infection()
            self.save_network()
            while self.simtime < self.exittime:
                self.time_step()
                self.save_network()
                self.store_states()
                
                if not self.contagion_network.is_active():
                    self.start_infection(first_infection = False)
                
            self.exit_sim()
                    
        elif infection_mode == 'single_infection':
            self.start_infection()
            self.save_network()
            while self.simtime < self.exittime:
                self.time_step()
                self.save_network()
                self.store_states()
            
            self.exit_sim()

        elif infection_mode == 'constant_infection':
            self.start_infection()
            self.save_network()
            while self.simtime < self.exittime:
                self.time_step()
                self.save_network()
                self.store_states()
            
                if self.simtime % self.reinfectiontime == 0:
                    self.start_infection(first_infection = False)
            
            self.exit_sim()
                
                
            
            # include options for seeding new cascades. This could be:
            #        size of seeding event could change. (not supported)  
            # create data organization tool???