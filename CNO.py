# this file needs documentation

import xgi
import random
import pandas as pd
import numpy as np
import scipy as sp
import time
import matplotlib.pyplot as plt

class CNO:
    
    def __init__(self, H, opinions, infectionrate, recoveryrate, K, ep ,timescale = 0.25):
        self.network = H
        self.dt = timescale
        self.gamma = recoveryrate
        self.betabar = infectionrate
        self.K = K
        self.ep = ep
        nodes = self.network.nodes
        for i in nodes:
            self.network.add_node(i, opinion = opinions[i], state = 0)
        
    def seed_contagion(self,num_infect_p1, num_infect_n1, configuration = 'random', first_infection = True):    
        if configuration == 'random':
            if first_infection:
                nodelist = list(self.network.nodes)
                p1_nodes = random.sample(nodelist, num_infect_p1)
                for node in p1_nodes:
                    nodelist.remove(node)
                    #self.network.add_node(node, state = 1, opinion = self.K)
                    self.network.add_node(node, state = 1)
                
                n1_nodes = random.sample(nodelist, num_infect_n1)
                for node in n1_nodes:
                    #self.network.add_node(node, state = -1, opinion = -self.K)
                    self.network.add_node(node, state = -1)       

            else:
                p1_nodelist = list(self.network.nodes.filterby_attr('opinion',  0, mode='gt'))
                n1_nodelist = list(self.network.nodes.filterby_attr('opinion',  0, mode='lt'))
                if num_infect_p1 < len(p1_nodelist):
                    p1_nodes = random.sample(p1_nodelist, num_infect_p1)
                else:
                    p1_nodes = p1_nodelist

                if num_infect_n1 < len(n1_nodelist):
                    n1_nodes = random.sample(n1_nodelist, num_infect_n1)
                else:
                    n1_nodes = n1_nodelist

                for node in p1_nodes:
                    self.network.add_node(node, state = 1)
                
                for node in n1_nodes:
                    self.network.add_node(node, state = -1)

        else:
            for pair in configuration:
                self.network.add_node(pair[0], state = pair[1])
        
    def get_spreaders(self,num_spreaders):
        spreaders = random.sample(list(self.network.nodes), num_spreaders)
        return spreaders
    
    def spread_event(self,spreaders):
        statedict = self.network.nodes.attrs('state').asdict()
        opiniondict = self.network.nodes.attrs('opinion').asdict()
        for spreader_node in spreaders:
            spreader_state = statedict[spreader_node]
            if spreader_state != 0:

                exposed_node = random.sample(list(self.network.nodes.neighbors(spreader_node)), 1)[0]

                exposed_node_state = statedict[exposed_node]
                old_opinion = opiniondict[exposed_node]

                #delta = np.sign(spreader_state-old_opinion)
                delta = spreader_state-old_opinion
                new_opinion = old_opinion + self.K*delta*self.dt
                #if np.abs(new_opinion) > 1:
                #    new_opinion = spreader_state 

                self.network.add_node(exposed_node, opinion = new_opinion)

                if exposed_node_state == 0:
                    self.infection_check(exposed_node, statedict[exposed_node], new_opinion, spreader_state)
            
    def heal_event(self, infected):
        opiniondict = self.network.nodes.attrs('opinion').asdict()
        statedict = self.network.nodes.attrs('state').asdict()
        for node in infected:
            self.heal_check(node, statedict[node], opiniondict[node])
            
    def infection_check(self,exposed_node, exposed_node_state, exposed_node_opinion, spreader_state):
        if random.random() < self.get_betabar(exposed_node, spreader_state, exposed_node_opinion)*self.dt:
            self.network.add_node(exposed_node, state = spreader_state)
            
    def heal_check(self,node,state,opinion):
        if random.random() < self.get_gamma(node,state,opinion)*self.dt:
            self.network.add_node(node, state = 0)
    
    # not in use
    def get_suscep_neighbors(self, node):
        neighbors  = list(self.network.nodes.neighbors(node))
        statedict = self.network.nodes.attrs('state').asdict()
        for Nnode in neighbors:
            if statedict[Nnode] != 0:
                neighbors.remove(Nnode)        
        return neighbors
    
    def is_active(self):
        nodes = self.network.nodes
        statedict = nodes.attrs('state').asdict()
        state = False
        for node in nodes:
            if statedict[node]!=0:
                state = True # some nodes are still infected        
        return state
    
    def get_gamma(self,node,state,opinion):
        return ((np.abs(opinion-state)+self.ep)/(2+self.ep))*self.gamma
        #if state == 1:
        #    return ((-(opinion-state)+self.ep)/(2+self.ep))*self.gamma # consider a different function for gamma(node) and betabar(node). There is a sigularity that exists in the mean field with this choice.
        #else:
        #    return (((opinion-state)+self.ep)/(2+self.ep))*self.gamma

    def get_betabar(self,node,state,opinion):
        return ((np.abs(opinion+state)+self.ep)/(2+self.ep))*self.betabar
        #if state == 1:
        #    return (( (opinion+state)+self.ep)/(2+self.ep))*self.betabar
        #else:
        #    return (( -(opinion+state)+self.ep)/(2+self.ep))*self.betabar

    def get_opinion_avg(self,state = 'all'):
        if state == 'all':
            return sum(self.network.nodes.attrs('opinion').aslist())/self.network.num_nodes
        else:
            return sum(self.network.nodes.filterby_attr('state', state).attrs('opinion').aslist())/self.network.num_nodes

    def get_frac_state(self,state):
        return len(self.network.nodes.filterby_attr('state', state))/self.network.num_nodes