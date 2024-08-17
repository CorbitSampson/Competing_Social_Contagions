import numpy as np
import math as math
import matplotlib.pyplot as pyplot

class meanfield_object:

    def __init__(self, infectionrate, recoverrate, C, ep, m, meandeg, initial_condition):
        self.gamma = recoverrate
        self.betabar = infectionrate
        self.K = C
        self.ep = ep
        self.m = m
        self.mk = meandeg
        self.v0 = initial_condition

    def odefunc(self, t, v):
        vprime = np.zeros(3)
        vprime[0] = -v[0]*self.get_gamma(v[2],1)  + self.m*(self.mk-1)/self.mk*(1 - v[0] - v[1])*v[0]*self.get_betabar(v[2],1)
        vprime[1] = -v[1]*self.get_gamma(v[2],-1) + self.m*(self.mk-1)/self.mk*(1 - v[0] - v[1])*v[1]*self.get_betabar(v[2],-1)
        vprime[2] = self.K*self.m*(np.sign(1-v[2])*v[0] + np.sign(-1 - v[2])*v[1])
        #vprime[2] = self.K*self.m*((1-v[2])*v[0] + (-1 - v[2])*v[1])


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


    def get_gamma(self,x,s):
        return self.gamma*(np.abs(x-s)+self.ep)/(2+self.ep)
        #if s == 1:
        #    return self.gamma*((s-x)+self.ep)/(2+self.ep)
        #else:
        #    return self.gamma*((x-s)+self.ep)/(2+self.ep)

    def get_betabar(self,x,s):
        return self.betabar*(np.abs(x+s)+self.ep)/(2+self.ep)
        #if s == 1:
        #    return self.betabar*((x+s)+self.ep)/(2+self.ep)
        #else:
        #    return self.betabar*(-(x+s)+self.ep)/(2+self.ep)

    def euler_solver(self, tspan, h):
        sys_size = len(self.v0)
        num_steps = math.ceil((tspan[1] - tspan[0])/h)
        t = [tspan[0] + i*h for i in range(0,num_steps)]
        y = np.zeros((sys_size, len(t)))
        y[:,0] = self.v0
        for i in range(0,num_steps-1):
            y[:,i+1] = y[:,i] + h*self.odefunc(t[i], y[:,i])
        return t,y

    
