import numpy as np
from numpy import random

def normal_wrapped(mu,sigma,n):
    
    samples =np.random.normal(mu,sigma,n)
    
    num_OL = len(samples[np.abs(samples) > 1])
    while num_OL > 0:
        samples[np.abs(samples) > 1] = np.random.normal(mu,sigma, num_OL)
        num_OL = len(samples[np.abs(samples) > 1])
        
    return samples

def delta_bimodal(mu,diff,n):

    val1 = mu + diff/2
    val2 = mu - diff/2
    
    samples = np.random.choice([val1,val2], n)
    return samples