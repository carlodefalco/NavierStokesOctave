#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: claudiosanavio
class of Grad-Carleman
"""

import numpy as np
import pandas as pd
import os,sys

from scipy.sparse import csr_matrix, hstack, vstack, eye #sparse matrices to save space

class grad_non_local:
    
    def __init__(self,L,dt,dx):

        self.length = L
        self.sites = L**2
        self.omega = 1
        self.cs = 1

        self.J = None #first order
        self.K = None #second ...
        self.L = None

        self.dt = dt
        self.dx = dx

        self.D = None
        
        self.order = 0

        self.A = None
        self.B = None
        self.C = None

    def initialize(self, Carleman = True, order = 2):
    
    #Here we define the first derivatives on the dim directions and the laplacian 
        Dx = (np.diag([0.5/self.dx]*(self.length-1),1)-np.diag([0.5/self.dx]*(self.length-1),-1))
        Dx[0,self.length-1]=-0.5/self.dx
        Dx[self.length-1,0]=0.5/self.dx

        self.D = [eye(self.sites)]
        for i in range(2):
            Di = 1
            for j in range(2):
                if(j==i):   Di = np.kron(Dx,Di)
                else:       Di = np.kron(np.eye(self.length),Di)
            self.D.append(csr_matrix(Di))

        self.order = order

        if(Carleman == True):     
            
            self.A = []
            for i in range(7):
                self.A.append([])
                for j in range(7):
                    self.A[-1].append([])
            self.A[0][0] = eye(self.sites)
            self.A[0][1] = self.dt*self.D[0]
            self.A[0][2] = self.dt*self.D[1]

            self.A[1][1] = eye(self.sites)
            self.A[1][3] = self.dt*self.D[0]
            self.A[1][4] = self.dt*self.D[1]

            self.A[2][2] = eye(self.sites)
            self.A[2][5] = self.dt*self.D[0]
            self.A[2][6] = self.dt*self.D[1]

            self.A[3][0] = self.omega*self.dt*self.vt*eye(self.sites)
            self.A[3][1] = 3*self.dt*self.D[0]
            self.A[3][2] = self.dt*self.D[1]
            self.A[3][3] = (1-self.dt*self.omega)*eye(self.sites)

            self.A[4][1] = self.dt*self.D[1]
            self.A[4][2] = self.dt*self.D[0]
            self.A[4][4] = (1-self.dt*self.omega)*eye(self.sites)

            self.A[5][1] = self.dt*self.D[1]
            self.A[5][2] = self.dt*self.D[0]
            self.A[5][5] = (1-self.dt*self.omega)*eye(self.sites)

            self.A[6][0] = self.omega*self.dt*self.vt*eye(self.sites)
            self.A[6][1] = self.dt*self.D[0]
            self.A[6][2] = 3*self.dt*self.D[1]
            self.A[6][6] = (1-self.dt*self.omega)*eye(self.sites)

            self.B = []
            for i in range(7):
                self.B.append([])
                for j in range(7):
                    self.B[-1].append([])                 
                    for k in range(7):
                        self.B[-1][-1].append([])            
            
            self.B[3][1][1] = 2*self.dt*self.omega*eye(self.sites)

            self.B[4][1][2] = 2*self.dt*self.omega*eye(self.sites)

            self.B[5][2][1] = 2*self.dt*self.omega*eye(self.sites)

            self.B[6][2][2] = 2*self.dt*self.omega*eye(self.sites)

            self.C = []
            for i in range(7):
                self.C.append([])
                for j in range(7):
                    self.C[-1].append([])                 
                    for k in range(7):
                        self.C[-1][-1].append([])  
                        for l in range(7):
                            self.C[-1][-1][-1].append([])  
                        
            self.C[3][1][1][1] = 4*self.dt*self.D[0]
            self.C[3][1][1][2] = 4*self.dt*self.D[1]
            
            self.C[4][1][2][1] = 4*self.dt*self.D[0]
            self.C[4][1][2][2] = 4*self.dt*self.D[1]

            self.C[5][2][1][1] = 4*self.dt*self.D[0]
            self.C[5][2][1][2] = 4*self.dt*self.D[1]

            self.C[6][2][2][1] = 4*self.dt*self.D[0]
            self.C[6][2][2][2] = 4*self.dt*self.D[1]
    
    def coord(self,x): #from the site index gives back its coordinates
        coord = list(np.base_repr(x, base=self.length, padding=3))[::-1]
        return coord[:3]

# initialization of kolmogorov-like flow
    def set_kolmogorov(self, amp = [0.1,0], k = [1,1], Carleman = True, NSE = True):  

            kappa = [2*np.pi/self.length*k[i] for i in range(len(k))] # kolmogorov flow

            self.J = []
            for q in range(7):  self.J.append([])
            self.J[0] = np.ones(self.sites)
            self.J[1] = []
            self.J[2] = []
            for x in range(self.sites): 
                self.J[1].append(amp[0]*np.cos(kappa[0]*int(self.coord(x)[0],self.length)))
                self.J[2].append(amp[1]*np.cos(kappa[1]*int(self.coord(x)[1],self.length)))
            self.J[1] = np.asarray(self.J[1])
            self.J[2] = np.asarray(self.J[2])

            if(NSE == True):
                self.J[3] = self.J[1]*self.J[1]+self.cs
                self.J[4] = self.J[1]*self.J[2]
                self.J[5] = self.J[1]*self.J[2]
                self.J[6] = self.J[2]*self.J[2]+self.cs
            else:
                self.J[3] = np.ones(self.sites)
                self.J[4] = np.zeros(self.sites)
                self.J[5] = np.zeros(self.sites)
                self.J[6] = np.ones(self.sites)

            if(Carleman == True):
                self.K = []
                self.L = []
                for i in range(7):
                    self.K.append([])
                    self.L.append([])
                    for j in range(7):
                        self.K[i].append(np.kron(self.J[i],self.J[j]))
                        self.L[-1].append([])
                        for k in range(7):
                            self.L[i][j].append(np.kron(np.kron(self.J[i],self.J[j]),self.J[k]))

    def set_Carleman(self, J_in = None, K_in = None, L_in = None, M_in = None, N_in = None):
            self.J = J_in
            self.K = K_in
            self.L = L_in

    def evolve(self, Carleman = True): 
        # Carleman evolution. new_J= ctoJ+JtoJ.J+KtoJ.K+LtoJ.L and analogously for the other Carleman variables
        # set_Carleman updates the variables
        if(Carleman == True):
            if(self.order == 1):
                new_J = []
                for i in range(7):
                    new_J.append(0)
                    for j in range(7):
                        new_J[j] += self.A[i][j]@self.J[j]
                self.set_Carleman(J_in = new_J)

            elif(self.order == 2):
                new_J = []
                new_K = []
                for i in range(7):
                    new_J.append(0)
                    new_K.append([])
                    for j in range(7):
                        new_J[i] += self.A[i][j]@self.J[j]
                        new_K[-1].append(0)
                        for k in range(7):
                            new_J[i] += self.B[i][j][k]@self.K[j][k]
                            for l in range(7):
                                new_K[i][j] += self.A[i][k]@self.A[j][l]@self.K[k][l]
                self.set_Carleman(J_in = new_J, K_in = new_K)

            elif(self.order == 3):
                new_J = []
                new_K = []
                new_L = []
                for i in range(7):
                    new_J.append(0)
                    new_K.append([])
                    new_L.append([])
                    for j in range(7):
                        new_J[i] += self.A[i][j]@self.J[j]
                        new_K[-1].append(0)
                        new_L[-1].append([])
                        for k in range(7):
                            new_J[i] += self.B[i][j][k]@self.K[j][k]
                            new_L[-1][-1].append(0)
                            for l in range(7):
                                new_J[i] += self.C[i][j][k]@self.L[j][k][l]
                                new_K[i][j] += self.A[i][k]@self.A[j][l]@self.K[k][l]
                                for m in range(7):
                                    new_K[i][j] += self.A[i][k]@self.B[j][l][m]@self.L[k][l][m]
                                    for n in range(7):
                                        new_L[i][j][k] += self.A[i][l]@self.A[j][m]@self.A[k][n]@self.L[l][m][n]
                self.set_Carleman(J_in = new_J, K_in = new_K, L_in = new_L)

        else: # Direct evolution
            new_J = self.J.copy()
            for i in range(3):
                new_J[i] += self.dt*(self.D[0].dot(self.J[2*i+1])+self.D[1].dot(self.J[2*i+2])) #rho, J1, J2
            new_J[3] += self.dt*(-self.omega*self.J[3]
                                 -self.omega*self.J[1]*self.J[1]*(2-self.J[0])
                                 -self.cs**2*self.omega*self.J[0]
                                -3*self.cs**2*self.D[0].dot(self.J[1])+self.D[1].dot(self.J[2]))
            new_J[4] += self.dt*(-self.omega*self.J[4]
                                 -self.omega*self.J[1]*self.J[2]*(2-self.J[0])
                                -self.cs**2*self.D[1].dot(self.J[1])
                                -self.cs**2*self.D[0].dot(self.J[2]))
            new_J[5] += self.dt*(-self.omega*self.J[5]
                                 -self.omega*self.J[1]*self.J[2]*(2-self.J[0])
                                -self.cs**2*self.D[1].dot(self.J[1])
                                -self.cs**2*self.D[0].dot(self.J[2]))
            new_J[6] += self.dt*(-self.omega*self.J[6]
                                 -self.omega*self.J[2]*self.J[2]*(2-self.J[0])
                                 -self.cs**2*self.omega*self.J[0]
                                -self.cs**2*self.D[0].dot(self.J[1])+3*self.D[1].dot(self.J[2]))
           
    def print_to_file(self, filepath):
        for x in range(self.sites):
            if(x == 0): df = pd.DataFrame([({**{f'x_{d}': int(self.coord(x)[d],self.length) for d in range(2)},**{f'J_{n}':self.J[n][x] for n in range(len(self.J))}})])
            else: df = pd.concat([df,pd.DataFrame([({**{f'x_{d}': int(self.coord(x)[d],self.length) for d in range(2)},**{f'J_{n}':self.J[n][x] for n in range(len(self.J))}})])])
        df.to_csv(filepath, float_format='%.8f',mode= 'a', index=False, header=not os.path.exists(filepath))






