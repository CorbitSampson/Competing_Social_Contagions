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