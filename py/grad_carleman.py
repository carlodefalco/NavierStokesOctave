#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: claudiosanavio

NS kolmogorov-like flow simulated using Carleman
"""

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import pandas as pd
import grad_class

#Nx,Ny,Nz = 32,32,1
L = 32 # sidelength of the lattice
dt,dx = 0.0001,0.1 #time and space discretization

#initial amplitudes of Kolmogorov flow
amp_x = 0.1 
amp_y = 0.

# wave-numbers of kolmogorov flow
k_x = 1
k_y = 1

#simulation time, always use ti=0 when Carleman = True
ti = 0
tf = 10000

# Carleman or DNS
Carleman = False
order = 1

#create class with relevant parameters
model = grad_class.grad_non_local(L,dt,dx)

#initialize variables
model.initialize(Carleman = Carleman, order = order)

if(Carleman == True):   dirname = 'Carleman_'+str(Carleman)+'_'+str(int(order))+'/L_'+str(L)+'_dx_'+str(int(dx*10000))+'_dt_'+str(int(dt*10000))+'_ax_'+str(int(amp_x*100))+'_ay_'+str(int(amp_y*100))+'_kx_'+str(int(k_x))+'_ky_'+str(int(k_y))
else: dirname = 'DNS/L_'+str(L)+'_dx_'+str(int(dx*10000))+'_dt_'+str(int(dt*10000))+'_ax_'+str(int(amp_x*100))+'_ay_'+str(int(amp_y*100))+'_kx_'+str(int(k_x))+'_ky_'+str(int(k_y))
if not os.path.exists(dirname):    os.makedirs(dirname)

if (ti != 0): 
    J_in =  np.asarray(pd.read_csv(dirname+'/flow_t_'+str(ti)+'.csv')[[f'J_{n}' for n in range(model.dim)]])
    model.set_kolmogorov(t_in = ti, J_in=J_in) #initialize kolmogorov flow
else: 
    model.set_kolmogorov(amp = [amp_x,amp_y] , k = [k_x,k_y], Carleman=Carleman)
    model.print_to_file(dirname+'/flow_t_'+str(ti)+'.csv')

for t in range(ti+1,tf):
    model.evolve(Carleman = Carleman) #evolution (Carleman or DNS)
    model.print_to_file(dirname+'/flow_t_'+str(t)+'.csv')
