# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:42:03 2017

@author: noor0021
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 22 10:08:57 2017

@author: mcca0206
"""
import gzip
import numpy as np
import scipy as sp
import pylab as plt
from HotRod_class import *
from DREAM import *
import csv
   
    
# Read hotrod time as a class. note: ./ needed for path directory
fn_name = "../inputs/templog_H_R2.csv"
sensor_fn = "../inputs/hotrod_sensor_array_updated.dat"
param_fn = "../inputs/param_data_opt2.dat"
htrd = HotRod(fn_name, sensor_fn,param_fn, paramChoice = 2, radius_inner = 0.028)

# set up markov chain; play around with burn in and chains to get better fit
nchains = 10
npars = htrd.npar 
D = DREAM(nchains, nburn = 300,npairs = 1)
D.sampler(htrd)
    