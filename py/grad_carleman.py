#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: claudiosanavio

Grad's moment formulation of hydrodynamics 
simulated using Euler-Carleman
"""

# import numpy as np
import os
# import time
# import matplotlib.pyplot as plt
import pandas as pd
import grad_class

#Nx,Ny,Nz = 32,32,1
L = 32 # sidelength of the lattice
dt,dx = 1,1 #time and space discretization

omega = 0.4
#initial amplitudes of Kolmogorov flow
amp_x = 0.1
amp_y = 0.1

# wave-numbers of kolmogorov flow
k_x = 1
k_y = 1

#simulation time, always use ti=0 when Carleman = True
ti = 0
tf = 20

# Carleman or DNS
locality = True

#create class with relevant parameters
model = grad_class.grad_carleman(L,dt,dx,omega,locality=locality)


#initialize variables 
for order in [0]: #Carleman evolution defined up to order 5. If order = 0, the DNS with Euler is performed.
        Carleman = bool(order)
        model.initialize(Carleman = Carleman, order = order)

        if(Carleman == True):   dirname = 'Carleman_'+str(Carleman)+'_'+str(int(order))+'/L_'+str(L)+'_dx_'+str(int(dx*10000))+'_dt_'+str(int(dt*10000))+'_ax_'+str(int(amp_x*100))+'_ay_'+str(int(amp_y*100))+'_kx_'+str(int(k_x))+'_ky_'+str(int(k_y))+'_omega_'+str(int(100*omega))
        else: dirname = 'DNS/L_'+str(L)+'_dx_'+str(int(dx*10000))+'_dt_'+str(int(dt*10000))+'_ax_'+str(int(amp_x*100))+'_ay_'+str(int(amp_y*100))+'_kx_'+str(int(k_x))+'_ky_'+str(int(k_y))+'_omega_'+str(int(100*omega))
        if not os.path.exists(dirname):    os.makedirs(dirname)

        model.set_kolmogorov(amp = [amp_x,amp_y] , k = [k_x,k_y], Carleman=Carleman, NSE = True) #initialize Kolmogorov-like flow, NSE = True when P_ab(t=0) is the value at equilibrium
        model.print_to_file(dirname+'/flow_t_'+str(ti)+'.csv')

        for t in range(ti+1,tf):
            model.evolve(Carleman = Carleman) #evolution (Carleman or DNS)
            model.print_to_file(dirname+'/flow_t_'+str(t)+'.csv')
