import numpy as np
import math as math
import matplotlib.pyplot as pyplot

class reduced_meanfield_object:

    def __init__(self, r0, K, ep, timeconstant, initial_condition, forcing_type = 'none', forcing_period = 0, forcing_amp = 0):
        self.r0 = r0
        self.K = K
        self.ep = ep
        self.TC = timeconstant
        self.v0 = initial_condition
        self.forcing_type = forcing_type
        self.forcing_period = forcing_period
        self.forcing_amp = forcing_amp

    def odefunc(self, t, v):
        vprime = np.zeros(3)
        vprime[0] = -v[0]*(np.abs(v[2]-1)+self.ep) + self.r0*(1-v[0]-v[1])*(v[0])*(np.abs(v[2]+1)+self.ep)
        vprime[1] = -v[1]*(np.abs(v[2]+1)+self.ep) + self.r0*(1-v[0]-v[1])*(v[1])*(np.abs(v[2]-1)+self.ep)
        #if np.abs(v[2]) < 1:
        #    vprime[2] = self.K*(np.sign(1-v[2])*v[0] + np.sign(-1 - v[2])*v[1])
        #else:
        #    vprime[2] = 0
        vprime[2] = self.K*((1-v[2])*v[0] + (-1 - v[2])*v[1])

        #vprime[0] = -v[0]*self.get_gamma(v[3],1)  + (self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[0]*self.get_betabar(v[2],1)
        #vprime[1] = -v[1]*self.get_gamma(v[4],-1) + (self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[1]*self.get_betabar(v[2],-1)
        #vprime[2] = (self.mk-1)/self.mk*self.K*self.m*(1-v[0] - v[1])*((1-v[2])*v[0] + (-1 - v[2])*v[1]) \
        #          + v[3]*(v[0]*self.get_gamma(v[2],1)) + v[4]*(v[1]*self.get_gamma(v[2],-1)) \
        #          - v[3]*((self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[0]*self.get_betabar(v[2],1)) \
        #          - v[4]*((self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[1]*self.get_betabar(v[2],-1))
        #vprime[3] = (self.mk-1)/self.mk*self.K*self.m*(v[0])*((1-v[3])*v[0] + (-1 - v[3])*v[1]) \
        #          - v[3]*(v[0]*self.get_gamma(v[2],1)) + v[3]*((self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[0]*self.get_betabar(v[2],1))
        #vprime[4] = (self.mk-1)/self.mk*self.K*self.m*(v[1])*((1-v[4])*v[0] + (-1 - v[4])*v[1]) \
        #          - v[4]*(v[1]*self.get_gamma(v[2],-1)) + v[4]*((self.mk-1)/self.mk*self.m*(1 - v[0] - v[1])*v[1]*self.get_betabar(v[2],-1))
        return vprime
    
    def evalfunc(self,t,v):
        vprime = np.zeros(v.shape)
        for i in range(0,v.shape[1]-1):
            vprime[:,i] = self.odefunc(t[i], v[:,i])
        return vprime

    def reintro_forcing(self,tvec,h):
        sys_size = len(self.v0)
        forcing_vec = np.zeros((sys_size,len(tvec)))
        num_multi = int(tvec[-1]//self.forcing_period)
        for i in range(0,len(tvec)):
            force = 0
            for j in range(1,num_multi+1): 
                force = force + (self.forcing_amp*self.forcing_period)*self.ddelta(j*self.forcing_period,tvec[i],h)
            forcing_vec[1,i] = force

        return forcing_vec
    
    def constant_forcing(self, value):
        return value

    def ddelta(self,shift,t,h):
        if (t < shift + h/2) and (t > shift - h/2): 
            return 1
        else:
            return 0
    
    def euler_solver(self, tspan, h, value = 0):
        sys_size = len(self.v0)
        num_steps = math.ceil((tspan[1] - tspan[0])/h)
        t = [tspan[0] + i*h for i in range(0,num_steps)]
        y = np.zeros((sys_size, len(t)))
        y[:,0] = self.v0
        if self.forcing_type == 'none':
            for i in range(0,num_steps-1):
                y[:,i+1] = y[:,i] + h*self.odefunc(t[i], y[:,i])
            return t,y

        elif self.forcing_type == 'single_periodic':
            forcing_vec = self.reintro_forcing(t,h)
            for i in range(0,num_steps-1):
                y[:,i+1] = y[:,i] + h*(self.odefunc(t[i], y[:,i])) + (1 - y[0,i] - y[1,i] )*forcing_vec[:,i]
            return t,y

        elif self.forcing_type == 'constant_function':
            forcing_vec = np.zeros((sys_size))
            forcing_vec[1] = self.constant_forcing(value)
            for i in range(0,num_steps-1):
                y[:,i+1] = y[:,i] + h*(self.odefunc(t[i], y[:,i]) + (1 - y[0,i] - y[1,i] )*forcing_vec)
            return t,y

    
