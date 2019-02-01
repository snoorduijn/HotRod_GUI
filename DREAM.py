# -*- coding: utf-8 -*-
"""
Created on Tue May  9 20:32:52 2017

@author: James
"""

import numpy as np
import scipy as sp
from HotRod_class_v10 import *
import copy

def likelihood(model,chain):
    mult = 1
    err = model.temp_error*mult
    L = 0.    
    obs = model.temp_obs*mult #obs
    sim = model.temp_sim*mult #sim
    LL = -1./(2.*err**2) * ((obs-sim)**2.) # log likelihood
    L = sum(np.ravel(LL))
    for n in range(chain.npars):
        if chain.uniform[n] == True:
            if chain.proposal[n] < (chain.mean[n] - chain.width[n]) or chain.proposal[n] > (chain.mean[n] + chain.width[n]):
                L += -9999
                print(n, chain.proposal[n], chain.mean[n]- chain.width[n], chain.mean[n] + chain.width[n],'here')
        else:
            L += -1./(2.*chain.width[n]**2) * ((chain.mean[n]-chain.proposal[n])**2.)
    L = np.asscalar(L)
    return(L)
    
def logLikelihood(model,chain):
    err = np.ravel(model.temp_error)
    obs = np.ravel(model.temp_obs) #obs
    sim = np.ravel(model.temp_sim) #sim
    LL = -(err.__len__() / 2) * np.log(2 * np.pi) - (1 / 2) * np.nansum(np.log(err**2)) - (err.__len__() / (2 * err**2)) * np.sum(
        (obs - sim) ** 2) # log likelihood
    LL = np.sum(LL)
    for n in range(chain.npars):
        if chain.uniform[n] == True:
            if chain.proposal[n] < (chain.mean[n] - chain.width[n]) or chain.proposal[n] > (chain.mean[n] + chain.width[n]):
                LL += -9999999
        else:
            LL +=  np.sum(-(err.__len__() / 2) * np.log(2 * np.pi) - (err.__len__() / 2) * np.nansum(np.log(chain.width[n])) - (1
                    / (2 * 2.*chain.width[n])) * np.sum((chain.mean[n]-chain.proposal[n]) ** 2))
    LL = np.asscalar(LL)
    return(LL)


class DREAM:
    def __init__(self, nchains, nburn = 100, npairs = 1,prior = False):
        self.nc = nchains
#        self.npars = npars
        self.prior = prior
        self.CR = 0.5
        self.chains = []
        self.genr = 0
        self.delta = npairs
        self.nb = nburn
        self.burn = True
        self.ct = 0
        
    def create_chains(self):
        for i in range(self.nc):
            self.chains.append(chain(self.npars,prior = self.prior))
            
    def par_set(self,log,uniform,mean,width,pmin,pmax):
        for i in range(self.nc):
            self.chains[i].par_set(log,uniform,mean,width,pmin,pmax,rstart = True)
            
    def gen_mod(self):
        self.genr += 1
        if self.genr == 5:
            self.genr = 0
            
    def par_reset(self):
        for i in range(self.nc):
            self.chains[i].pars = self.chains[i].current       
            
    def propgen(self):
        if self.burn == True:
            self.dumloc = np.random.randint(self.ncr) 
            self.CR = (self.dumloc +1) / self.ncr
        for i in range(self.nc):
            #Identify random chains to use
            dum = list(range(self.nc))
            dum.remove(i)
            dx = np.zeros(self.npars)
            ll = self.nc-1
            for k in range(self.delta):
                dum2  = np.random.randint(ll)
                ll = ll-1
                dx += self.chains[dum[dum2]].current
                #print(dum2,dx)
                dum.remove(dum[dum2])
                dum2  = np.random.randint(ll)
                ll = ll-1
                dx += -self.chains[dum[dum2]].current
                #print(dum2,dx)
                dum.remove(dum[dum2])
            #print(dx)
            d = self.npars
            #print(d)
            mult = np.ones(self.npars)
            for k in range(self.npars):
                if sp.rand() > self.CR:
                    dx[k] = 0.
                    d += -1
                    mult[k] = 0.
            if self.genr < 4:
                if d == 0:
                    gamma = 0
                else:
                    gamma = 2.38/np.sqrt(2.*self.delta*d)
                if gamma > 1.0:
                    gamma = 1.0#                
                
            else:
                gamma = 1.0
            
            dx = dx * gamma # + np.random.normal(scale = self.chains[i].pwidth) * mult
            self.chains[i].proposal = self.chains[i].current + dx
            for k in range(self.npars):
                if self.chains[i].proposal[k] < self.chains[i].pmin[k]:
                     self.chains[i].proposal[k] = self.chains[i].pmax[k] - (self.chains[i].pmin[k] - self.chains[i].proposal[k])
                elif self.chains[i].proposal[k] > self.chains[i].pmax[k]:
                     self.chains[i].proposal[k] = self.chains[i].pmin[k] - (self.chains[i].pmax[k] - self.chains[i].proposal[k])
        self.ct+=1
              
    def delm_update(self):
        for i in range(self.nc):
            for j in range(self.npars):
                self.delm[self.dumloc] += (self.chains[i].pars[-1,j] - self.chains[i].pars[-2,j]) **2. /self.Var[j]                
                self.crct[self.dumloc] +=1
                
    def Chain_removal(self,):
        self.reset = False
        Omega = np.zeros((self.nc))
        n = int(self.ct/2)
        last = np.zeros((self.nc))
        for i in range(self.nc):
            Omega[i] = np.average(self.chains[i].likelihood[n:])
            last[i] = self.chains[i].likelihood[-1]
        best = np.argmax(last)
        IQRmin = np.percentile(Omega,25) - 2 * (np.percentile(Omega,75) - np.percentile(Omega,25))        
        for i in range(self.nc):      
            if Omega[i] < IQRmin:
                print("Chain %i replaced (%.2f => %.2f)"%(i, Omega[i],IQRmin))
                self.chains[i].current = np.copy(self.chains[best].current)
                self.chains[i].Lold = np.copy(self.chains[best].Lold)
                self.reset = True
        if self.reset == False and self.ct > self.nb:
            self.burn = False
        elif self.reset == True:
            self.reset = False
            self.ct = -1
            self.genr = 0
            for i in range(self.nc):
                self.chains[i].likelihood = [self.chains[i].Lold]
                self.chains[i].pars = self.chains[i].pars[-1,:]

    def set_CR(self, ncr = 3):
        self.crct=np.zeros(ncr)
        self.ncr = ncr
        self.delm = np.zeros(ncr)
        
    def Rget(self):
        mean = np.zeros((self.nc,self.npars))
        var = np.zeros((self.nc,self.npars))
        for i in range(self.nc):
            mean[i,:], var[i,:], n = self.chains[i].varget() 
        W = np.average(var,axis=0)
        W[W == 0.] = np.nan
        B = np.var(mean,axis=0) * n
        self.Var = (1. -1./n) * W + 1. /n * B
        self.R = np.sqrt(self.Var/W)
        
    def sampler(self, parent, convergence_criteria = 1.2):
        print(parent.paramChoice)
        self.npars = len(parent.Log)
        print(self.npars)
        self.create_chains()
        fn_name = parent.temp_fn
        
        self.par_set(parent.Log, parent.Unif, 
          parent.Mean, parent.width, 
          parent.pmin, parent.pmax)
        
        #make a list of models
        models = []
        for i in range(self.nc):
            child = copy.deepcopy(parent)
            child.ID = i
            models.append(child)

        for i in range(self.nc):
            models[i].set_propositions(self.chains[i].current)
            self.chains[i].proposal = np.copy(self.chains[i].current)
            models[i].get_set_solve()
            self.chains[i].Lold = likelihood(models[i],self.chains[i])
#            print(self.chains[i].Lold)    
        
        self.set_CR()
        
        while self.burn == True:
            self.propgen()
            
            for i in range(self.nc):
                models[i].set_propositions(self.chains[i].proposal)
                models[i].get_set_solve()
                self.chains[i].Lnew = likelihood(models[i], self.chains[i])
                self.chains[i].tprob() 
                
            if self.ct >= 2:
                for i in range(self.nc):
                    self.chains[i].varget()
                self.Rget()
                if self.ct%10 == 0:
                    print("Burn #%i, R = %.4f"%(self.ct, max(self.R)))
                self.delm_update()
            if self.ct >= 5:
                self.Chain_removal()
            self.gen_mod()
                     
        self.delm = self.delm/self.crct
        self.CR = (np.argmax(self.delm)+1)/self.ncr
        self.ct = 0
        self.genr = 0
        self.R[:] = 10      

        f = open('%s_param_data.dat'%(fn_name[:-4]), 'w') 
        g = open('%s_sim_data.dat'%(fn_name[:-4]), 'w') 
        for i in range(self.nc):
            self.chains[i].pars=self.chains[i].current
            self.chains[i].likelihood = [self.chains[i].Lold]
        
        while max(self.R) > convergence_criteria:
            self.propgen()
        
            for i in range(self.nc):
                models[i].set_propositions(self.chains[i].proposal)
                models[i].get_set_solve()
                self.chains[i].Lnew = likelihood(models[i],self.chains[i])
                self.chains[i].tprob() 
                self.chains[i].varget()
                
            if self.ct > 10:
                self.Rget()
            if self.ct%5 == 0:
                print("Sample %i, Max R stat: %.4f"%(self.ct, max(self.R)))
            self.gen_mod() 
        
        
            for i in range(self.nc):
                dim = np.shape(self.chains[i].pars)   
                f.write('%.15e ' % self.chains[i].likelihood[-1])
                for k in range(dim[1]):
                    f.write('%.7e ' % self.chains[i].pars[-1,k])
                f.write('\n') 
        print("Sampling Complete (total # samples %i): "%self.ct + " ".join(map(str,self.R)))
        f.close()
        g.close()
        
class chain:
    def __init__(self,npars,prior = False):
        self.npars = npars
        self.current = np.zeros(npars)
        self.proposal = np.zeros(npars)
        self.Lold = 0.
        self.Lnew = 0.
        self.pwidth = 0.
        self.prior = prior
        self.likelihood = []
        
    def par_set(self,log,uniform,mean,width,pmin,pmax,rstart = False):
        self.log = log #logical
        self.uniform = uniform # logical
        self.mean = mean # mean or centre of distribution (log for log parameters)
        self.width = width #stabdard deviation or width (;og for log)
        self.rs = rstart #logica - if you want random generte start point
        self.pwidth = []    
        self.pmin=pmin
        self.pmax=pmax
        for i in range(np.size(self.width)):        
            self.pwidth.append(self.width[i]/25.)
        if self.rs == True:
            for i in range(self.npars):
                if uniform[i] == True:
                    self.current[i] = self.mean[i] - self.width[i] + sp.rand() * 2. *self.width[i]
                else:
                    self.current[i] = np.random.normal(loc = self.mean[i], scale = self.width[i])
        else:
            self.current[i] = np.copy(self.mean[i])           
        self.pars = np.copy(self.current)
        
        
    def prior(self):
        for i in range(self.npars):
            if self.uniform == True:
                if self.proposal[i] < (self.mean[i] - self.width[i]) or self.proposal[i] > (self.mean[i] + self.width[i]):
                    self.Lnew += -1000
            else:
                self.Lnew += ((self.proposal[i]-self.mean[i])**2.)/(2.*self.width[i]**2.)
                
    def tprob(self):
        if self.Lnew > self.Lold:
            self.current = np.copy(self.proposal)
            self.Lold = np.copy (self.Lnew)
        else:
            if sp.rand() < np.exp(self.Lnew - self.Lold):
                self.current = np.copy(self.proposal)
                self.Lold = np.copy(self.Lnew)
        self.pars = np.vstack((self.pars,self.current))
        self.likelihood.append(self.Lold)
        
    def varget(self):
        n = int(np.shape(self.pars)[0]/2)
        s2 = np.var(self.pars[n:,:],axis = 0)
        mn = np.average(self.pars[n:,:], axis=0)
        return(mn,s2,n)
    
 