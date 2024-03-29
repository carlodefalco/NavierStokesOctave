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

class grad_carleman:

    def __init__(self,L,dt,dx,omega,locality = True):

        self.length = L
        self.sites = L**2
        self.omega = omega
        self.cs = 1/np.sqrt(3)

        self.J = None #first order
        self.K = None #second ...
        self.L = None
        self.M = None
        self.N = None
        self.O = None

        self.dt = dt
        self.dx = dx

        self.D = None
        
        self.order = 2

        self.A = None
        self.B = None
        self.C = None

        self.locality = locality

    def initialize(self, Carleman = True, order = 2):
    #Here we define the first derivatives on the dim directions and the laplacian 
        Dx = (np.diag([0.5/self.dx]*(self.length-1),1)-np.diag([0.5/self.dx]*(self.length-1),-1))
        Dx[0,self.length-1]=-0.5/self.dx
        Dx[self.length-1,0]=0.5/self.dx

        self.D = []
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
            self.A[0][1] = -self.dt*self.D[0]
            self.A[0][2] = -self.dt*self.D[1]

            self.A[1][1] = eye(self.sites)
            self.A[1][3] = -self.dt*self.D[0]
            self.A[1][4] = -self.dt*self.D[1]

            self.A[2][2] = eye(self.sites)
            self.A[2][5] = -self.dt*self.D[0]
            self.A[2][6] = -self.dt*self.D[1]

            self.A[3][0] = -self.omega*self.dt*self.cs**2*eye(self.sites)
            self.A[3][1] = -3.*self.dt*self.cs**2*self.D[0]
            self.A[3][2] = -self.dt*self.cs**2*self.D[1]
            self.A[3][3] = (1-self.dt*self.omega)*eye(self.sites)

            self.A[4][1] = -self.dt*self.cs**2*self.D[1]
            self.A[4][2] = -self.dt*self.cs**2*self.D[0]
            self.A[4][4] = (1-self.dt*self.omega)*eye(self.sites)


            self.A[5][1] = -self.dt*self.cs**2*self.D[1]
            self.A[5][2] = -self.dt*self.cs**2*self.D[0]
            self.A[5][5] = (1-self.dt*self.omega)*eye(self.sites)

            self.A[6][0] = self.omega*self.dt*self.cs**2*eye(self.sites)
            self.A[6][1] = -self.dt*self.cs**2*self.D[0]
            self.A[6][2] = -3.*self.dt*self.cs**2*self.D[1]
            self.A[6][6] = (1.-self.dt*self.omega)*eye(self.sites)

            self.B = []
            for i in range(7):
                self.B.append([])
                for j in range(7):
                    self.B[-1].append([])                 
                    for k in range(7):
                        self.B[-1][-1].append([])            
            
            self.B[3][1][1] = 2.*self.dt*self.omega*eye(self.sites)

            self.B[4][1][2] = 2.*self.dt*self.omega*eye(self.sites)

            self.B[5][2][1] = 2.*self.dt*self.omega*eye(self.sites)

            self.B[6][2][2] = 2.*self.dt*self.omega*eye(self.sites)


            self.C = []
            for i in range(7):
                self.C.append([])
                for j in range(7):
                    self.C[-1].append([])                 
                    for k in range(7):
                        self.C[-1][-1].append([])  
                        for l in range(7):
                            self.C[-1][-1][-1].append([])  
                        
            self.C[3][1][1][0] = -self.dt*self.omega
            
            self.C[4][1][2][0] = -self.dt*self.omega

            self.C[5][2][1][0] = -self.dt*self.omega

            self.C[6][2][2][0] = -self.dt*self.omega
    
    def coord(self,x): #from the site index gives back its coordinates
        coord = list(np.base_repr(x, base=self.length, padding=3))[::-1]
        return coord[:3]

# initialization of kolmogorov-like flow
    def set_kolmogorov(self, amp = [0.1,0], k = [1,1], Carleman = True, NSE = False):  

            kappa = [2*np.pi/(self.length)*k[i] for i in range(len(k))] # kolmogorov flow

            self.J = []
            for q in range(7):  self.J.append([])
            self.J[0] = np.ones(self.sites)
            self.J[1] = []
            self.J[2] = []
            for x in range(self.sites): 
                self.J[1].append(amp[0]*np.cos(kappa[0]*int(self.coord(x)[1],self.length)))
                self.J[2].append(amp[1]*np.cos(kappa[1]*int(self.coord(x)[0],self.length)))
            self.J[1] = np.asarray(self.J[1])
            self.J[2] = np.asarray(self.J[2])

            if(NSE == True):
                self.J[3] = self.J[1]*self.J[1]+self.cs**2*self.J[0]
                self.J[4] = self.J[1]*self.J[2]
                self.J[5] = self.J[1]*self.J[2]
                self.J[6] = self.J[2]*self.J[2]+self.cs**2*self.J[0]
            else:
                self.J[3] = np.ones(self.sites)
                self.J[4] = np.zeros(self.sites)
                self.J[5] = np.zeros(self.sites)
                self.J[6] = np.ones(self.sites)

            if(Carleman == True):
                if(self.order>1):
                    self.K = []
                    self.L = []
                    self.M = []
                    self.N = []
                    self.O = []
                    for i in range(7):
                        self.K.append([])
                        self.L.append([])
                        self.M.append([])
                        self.N.append([])
                        self.O.append([])
                        for j in range(7):
                            self.K[i].append(self.J[i]*self.J[j])
                            self.L[-1].append([])
                            self.M[-1].append([])
                            self.N[-1].append([])
                            self.O[-1].append([])
                            if(self.order>2):
                                for k in range(7):
                                    self.L[i][j].append(self.J[i]*self.J[j]*self.J[k])
                                    self.M[-1][-1].append([])
                                    self.N[-1][-1].append([])
                                    self.O[-1][-1].append([])
                                    if(self.order>3):
                                        for l in range(7):
                                            self.M[i][j][k].append(self.J[i]*self.J[j]*self.J[k]*self.J[l])
                                            self.N[-1][-1][-1].append([])
                                            self.O[-1][-1][-1].append([])
                                            if(self.order>4):
                                                for m in range(7):
                                                    self.N[i][j][k][l].append(self.J[i]*self.J[j]*self.J[k]*self.J[l]*self.J[m])
                                                    self.O[-1][-1][-1][-1].append([])
                                                    if(self.order>5):
                                                        for n in range(7):
                                                            self.O[i][j][k][l][m].append(self.J[i]*self.J[j]*self.J[k]*self.J[l]*self.J[m]*self.J[n])

                    # else:
                    #     self.K = []
                    #     # self.L = []
                    #     for i in range(7):
                    #         self.K.append([])
                    #         for j in range(7):
                    #             self.K[i].append(np.tensordot(self.J[i],self.J[j],axes = 0))
                            # self.L[-1].append([])
                            # for k in range(7):
                                # self.L[i][j].append(np.kron(np.kron(self.J[i],self.J[j]),self.J[k]))

    def set_Carleman(self, J_in = None, K_in = None, L_in = None, M_in = None, N_in = None, O_in = None):
            self.J = J_in
            self.K = K_in
            self.L = L_in
            self.M = M_in
            self.N = N_in
            self.O = O_in

    def evolve(self, Carleman = True): 

        if(Carleman == True):
            if(self.order == 1):
                new_J = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += (self.A[i][j]).dot(self.J[j])
                self.set_Carleman(J_in = new_J)

            elif(self.order == 2):
                new_J = []
                new_K = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    new_K.append([])
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += self.A[i][j].dot(self.J[j])
                        new_K[-1].append(np.zeros(self.sites))
                        for k in range(7):
                            if(self.B[i][j][k]!=[]):
                                new_J[i] += self.B[i][j][k].dot(self.K[j][k])
                            for l in range(7):
                                if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                                    new_K[i][j] += self.A[i][k].dot(self.A[j][l].dot(self.K[k][l]))
                self.set_Carleman(J_in = new_J, K_in = new_K)

            elif(self.order == 3):
                new_J = []
                new_K = []
                new_L = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    new_K.append([])
                    new_L.append([])
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += self.A[i][j].dot(self.J[j])
                        new_K[-1].append(np.zeros(self.sites))
                        new_L[-1].append([])
                        for k in range(7):
                            if(self.B[i][j][k]!=[]):
                                new_J[i] += self.B[i][j][k].dot(self.K[j][k])
                            for l in range(7):
                                if(self.C[i][j][k][l]!=[]):
                                    new_J[i]+= self.C[i][j][k][l]*self.L[j][k][l]
                                if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                                    new_K[i][j] += self.A[i][k].dot(self.A[j][l].dot(self.K[k][l])) 
                                new_L[-1][-1].append(np.zeros(self.sites))
                                for m in range(7):
                                    if(self.A[i][k]!=[] and self.B[j][l][m]!=[]):
                                        new_K[i][j] += self.A[i][k].dot(self.B[j][l][m].dot(self.L[k][l][m]))
                                    if(self.A[j][m]!=[] and self.B[i][k][l]!=[]):
                                        new_K[i][j] += self.B[i][k][l].dot(self.A[j][m].dot(self.L[k][l][m]))
                                    for n in range(7):
                                        if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.A[k][n]!=[]):
                                            new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.A[k][n].dot(self.L[l][m][n])))
                self.set_Carleman(J_in = new_J, K_in = new_K, L_in= new_L)

            elif(self.order == 4):
                new_J = []
                new_K = []
                new_L = []
                new_M = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    new_K.append([])
                    new_L.append([])
                    new_M.append([])
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += self.A[i][j].dot(self.J[j])
                        new_K[-1].append(np.zeros(self.sites))
                        new_L[-1].append([])
                        new_M[-1].append([])
                        for k in range(7):
                            if(self.B[i][j][k]!=[]):
                                new_J[i] += self.B[i][j][k].dot(self.K[j][k])
                            for l in range(7):
                                if(self.C[i][j][k][l]!=[]):
                                    new_J[i]+= self.C[i][j][k][l]*self.L[j][k][l]
                                if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                                    new_K[i][j] += self.A[i][k].dot(self.A[j][l].dot(self.K[k][l])) 
                                new_L[-1][-1].append(np.zeros(self.sites))
                                new_M[-1][-1].append([])
                                for m in range(7):
                                    if(self.A[i][k]!=[] and self.B[j][l][m]!=[]):
                                        new_K[i][j] += self.A[i][k].dot(self.B[j][l][m].dot(self.L[k][l][m]))
                                    if(self.A[j][m]!=[] and self.B[i][k][l]!=[]):
                                        new_K[i][j] += self.B[i][k][l].dot(self.A[j][m].dot(self.L[k][l][m]))
                                    new_M[-1][-1][-1].append(np.zeros(self.sites))
                                    for n in range(7):
                                        if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.A[k][n]!=[]):
                                            new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.A[k][n].dot(self.L[l][m][n])))
                                        if(self.A[i][k]!=[] and self.C[j][l][m][n]!=[]):
                                            new_K[i][j] += self.A[i][k].dot(self.C[j][l][m][n]*self.M[k][l][m][n])                        
                                        if(self.C[i][k][m][n]!=[] and self.A[j][l]!=[]):
                                            new_K[i][j] += self.C[i][k][m][n]*self.A[j][l].dot(self.M[k][l][m][n])                                      
                                        if(self.B[i][k][l]!=[] and self.B[j][m][n]!=[]):
                                            new_K[i][j] += self.B[i][k][l].dot(self.B[j][m][n].dot(self.M[k][l][m][n]))
                                        for o in range(7):
                                            if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.B[k][n][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.B[k][n][o].dot(self.M[l][m][n][o])))
                                            if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            for p in range(7):
                                                if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[] and self.A[l][p]!=[]):
                                                    new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.A[l][p].dot(self.M[m][n][o][p]))))
                self.set_Carleman(J_in = new_J, K_in = new_K, L_in= new_L, M_in= new_M)

            elif(self.order == 5):
                new_J = []
                new_K = []
                new_L = []
                new_M = []
                new_N = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    new_K.append([])
                    new_L.append([])
                    new_M.append([])
                    new_N.append([])
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += self.A[i][j].dot(self.J[j])
                        new_K[i].append(np.zeros(self.sites))
                        new_L[i].append([])
                        new_M[i].append([])
                        new_N[i].append([])
                        for k in range(7):
                            if(self.B[i][j][k]!=[]):
                                new_J[i] += self.B[i][j][k].dot(self.K[j][k])
                            new_L[i][j].append(np.zeros(self.sites))
                            new_M[i][j].append([])
                            new_N[i][j].append([])                            
                            for l in range(7):
                                if(self.C[i][j][k][l]!=[]):
                                    new_J[i]+= self.C[i][j][k][l]*self.L[j][k][l]
                                if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                                    new_K[i][j] += self.A[i][k].dot(self.A[j][l].dot(self.K[k][l])) 
                                new_M[i][j][k].append(np.zeros(self.sites))
                                new_N[i][j][k].append([])                                 
                                for m in range(7):
                                    if(self.A[i][k]!=[] and self.B[j][l][m]!=[]):
                                        new_K[i][j] += self.A[i][k].dot(self.B[j][l][m].dot(self.L[k][l][m]))
                                    if(self.A[j][m]!=[] and self.B[i][k][l]!=[]):
                                        new_K[i][j] += self.B[i][k][l].dot(self.A[j][m].dot(self.L[k][l][m]))
                                    new_N[i][j][k][l].append(np.zeros(self.sites))
                                    for n in range(7):
                                        if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.A[k][n]!=[]):
                                            new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.A[k][n].dot(self.L[l][m][n])))
                                        if(self.A[i][k]!=[] and self.C[j][l][m][n]!=[]):
                                            new_K[i][j] += self.A[i][k].dot(self.C[j][l][m][n]*self.M[k][l][m][n])                        
                                        if(self.C[i][k][m][n]!=[] and self.A[j][l]!=[]):
                                            new_K[i][j] += self.C[i][k][m][n]*self.A[j][l].dot(self.M[k][l][m][n])                                      
                                        if(self.B[i][k][l]!=[] and self.B[j][m][n]!=[]):
                                            new_K[i][j] += self.B[i][k][l].dot(self.B[j][m][n].dot(self.M[k][l][m][n]))
                                        for o in range(7):
                                            if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.B[k][n][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.B[k][n][o].dot(self.M[l][m][n][o])))
                                            if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            if(self.B[i][k][l]!=[] and self.C[j][m][n][o]!=[]):
                                                new_K[i][j] += self.C[j][m][n][o]*self.B[i][k][l].dot(self.N[k][l][m][n][o])
                                            if(self.C[i][k][l][m]!=[] and self.B[j][n][o]!=[]):
                                                new_K[i][j] += self.C[i][k][l][m]*self.B[j][n][o].dot(self.N[k][l][m][n][o])
                                            for p in range(7):
                                                if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[] and self.A[l][p]!=[]):
                                                    new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.A[l][p].dot(self.M[m][n][o][p]))))
                                                if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.B[k][o][p]!=[]):
                                                    new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.B[k][o][p].dot(self.N[l][m][n][o][p])))
                                                if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.B[k][o][p]!=[]):
                                                    new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.B[k][o][p].dot(self.N[l][m][n][o][p])))
                                                if(self.B[i][l][m]!=[] and self.B[j][n][o]!=[] and self.A[k][p]!=[]):
                                                    new_L[i][j][k] += self.B[i][l][m].dot(self.B[j][n][o].dot(self.A[k][p].dot(self.N[l][m][n][o][p])))
                                                for q in range(7):
                                                    if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[] and self.B[l][p][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.B[l][p][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.B[k][o][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.B[k][o][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.A[i][m]!=[] and self.B[j][n][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.B[j][n][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.B[i][m][n]!=[] and self.A[j][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.B[i][m][n].dot(self.A[j][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    for r in range(7):
                                                        if(self.A[i][n]!=[] and self.A[j][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[] and self.A[m][r]!=[]):
                                                            new_N[i][j][k][l][m] += self.A[i][n].dot(self.A[j][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.A[m][r].dot(self.N[n][o][p][q][r])))))                                                        
                self.set_Carleman(J_in = new_J, K_in = new_K, L_in= new_L, M_in= new_M, N_in= new_N)

            elif(self.order == 6):
                new_J = []
                new_K = []
                new_L = []
                new_M = []
                new_N = []
                new_O = []
                for i in range(7):
                    new_J.append(np.zeros(self.sites))
                    new_K.append([])
                    new_L.append([])
                    new_M.append([])
                    new_N.append([])
                    new_O.append([])
                    for j in range(7):
                        if(self.A[i][j]!=[]):
                            new_J[i] += self.A[i][j].dot(self.J[j])
                        new_K[i].append(np.zeros(self.sites))
                        new_L[i].append([])
                        new_M[i].append([])
                        new_N[i].append([])
                        new_O[i].append([])
                        for k in range(7):
                            if(self.B[i][j][k]!=[]):
                                new_J[i] += self.B[i][j][k].dot(self.K[j][k])
                            new_L[i][j].append(np.zeros(self.sites))
                            new_M[i][j].append([])
                            new_N[i][j].append([])                            
                            new_O[i][j].append([])                            
                            for l in range(7):
                                if(self.C[i][j][k][l]!=[]):
                                    new_J[i]+= self.C[i][j][k][l]*self.L[j][k][l]
                                if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                                    new_K[i][j] += self.A[i][k].dot(self.A[j][l].dot(self.K[k][l])) 
                                new_M[i][j][k].append(np.zeros(self.sites))
                                new_N[i][j][k].append([])                                 
                                new_O[i][j][k].append([])                                 
                                for m in range(7):
                                    if(self.A[i][k]!=[] and self.B[j][l][m]!=[]):
                                        new_K[i][j] += self.A[i][k].dot(self.B[j][l][m].dot(self.L[k][l][m]))
                                    if(self.A[j][m]!=[] and self.B[i][k][l]!=[]):
                                        new_K[i][j] += self.B[i][k][l].dot(self.A[j][m].dot(self.L[k][l][m]))
                                    new_N[i][j][k][l].append(np.zeros(self.sites))
                                    new_O[i][j][k][l].append([])
                                    for n in range(7):
                                        if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.A[k][n]!=[]):
                                            new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.A[k][n].dot(self.L[l][m][n])))
                                        if(self.A[i][k]!=[] and self.C[j][l][m][n]!=[]):
                                            new_K[i][j] += self.A[i][k].dot(self.C[j][l][m][n]*self.M[k][l][m][n])                        
                                        if(self.C[i][k][m][n]!=[] and self.A[j][l]!=[]):
                                            new_K[i][j] += self.C[i][k][m][n]*self.A[j][l].dot(self.M[k][l][m][n])                                      
                                        if(self.B[i][k][l]!=[] and self.B[j][m][n]!=[]):
                                            new_K[i][j] += self.B[i][k][l].dot(self.B[j][m][n].dot(self.M[k][l][m][n]))
                                        new_N[i][j][k][l][m].append(np.zeros(self.sites))
                                        for o in range(7):
                                            if(self.A[i][l]!=[] and self.A[j][m]!=[] and self.B[k][n][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.A[j][m].dot(self.B[k][n][o].dot(self.M[l][m][n][o])))
                                            if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[]):
                                                new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.M[l][m][n][o])))
                                            if(self.B[i][k][l]!=[] and self.C[j][m][n][o]!=[]):
                                                new_K[i][j] += self.C[j][m][n][o]*self.B[i][k][l].dot(self.N[k][l][m][n][o])
                                            if(self.C[i][k][l][m]!=[] and self.B[j][n][o]!=[]):
                                                new_K[i][j] += self.C[i][k][l][m]*self.B[j][n][o].dot(self.N[k][l][m][n][o])
                                            for p in range(7):
                                                if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[] and self.A[l][p]!=[]):
                                                    new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.A[l][p].dot(self.M[m][n][o][p]))))
                                                if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.B[k][o][p]!=[]):
                                                    new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.B[k][o][p].dot(self.N[l][m][n][o][p])))
                                                if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.B[k][o][p]!=[]):
                                                    new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.B[k][o][p].dot(self.N[l][m][n][o][p])))
                                                if(self.B[i][l][m]!=[] and self.B[j][n][o]!=[] and self.A[k][p]!=[]):
                                                    new_L[i][j][k] += self.B[i][l][m].dot(self.B[j][n][o].dot(self.A[k][p].dot(self.N[l][m][n][o][p])))
                                                if(self.C[i][k][l][m]!=[] and self.C[j][n][o][p]!=[]):
                                                    new_K[i][j] += self.C[i][k][l][m]*self.C[j][n][o][p]*self.O[k][l][m][n][o][p]
                                                for q in range(7):
                                                    if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.A[k][o]!=[] and self.B[l][p][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.A[k][o].dot(self.B[l][p][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.B[k][o][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.B[k][o][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.A[i][m]!=[] and self.B[j][n][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.A[i][m].dot(self.B[j][n][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.B[i][m][n]!=[] and self.A[j][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[]):
                                                        new_M[i][j][k][l] += self.B[i][m][n].dot(self.A[j][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.N[m][n][o][p][q]))))
                                                    if(self.A[i][l]!=[] and self.B[j][m][n]!=[] and self.C[k][o][p][q]!=[]):
                                                        new_L[i][j][k] += self.A[i][l].dot(self.B[j][m][n].dot(self.C[k][o][p][q]*self.O[l][m][n][o][p][q]))
                                                    if(self.B[i][l][m]!=[] and self.A[j][n]!=[] and self.C[k][o][p][q]!=[]):
                                                        new_L[i][j][k] += self.B[i][l][m].dot(self.A[j][n].dot(self.C[k][o][p][q]*self.O[l][m][n][o][p][q]))
                                                    if(self.B[i][l][m]!=[] and self.C[j][n][o][p]!=[] and self.A[k][q]!=[]):
                                                        new_L[i][j][k] += self.B[i][l][m].dot(self.C[j][n][o][p]*self.A[k][q].dot(self.O[l][m][n][o][p][q]))
                                                    
                                                    if(self.A[i][l]!=[] and self.C[j][m][n][o]!=[] and self.B[k][p][q]!=[]):
                                                        new_L[i][j][k] += self.A[i][l].dot(self.C[j][m][n][o]*self.B[k][p][q].dot(self.O[l][m][n][o][p][q]))                                                                                                                                                                    
                                                    if(self.C[i][l][m][n]!=[] and self.A[j][o]!=[] and self.B[k][p][q]!=[]):
                                                        new_L[i][j][k] += self.C[i][l][m][n]*self.A[j][o].dot(self.B[k][p][q].dot(self.O[l][m][n][o][p][q]))                                                                                                                                                                    
                                                    if(self.C[i][l][m][n]!=[] and self.B[j][o][p]!=[] and self.A[k][q]!=[]):
                                                        new_L[i][j][k] += self.C[i][l][m][n]*self.B[j][o][p].dot(self.A[k][q].dot(self.O[l][m][n][o][p][q]))                                                                                                                                                                   
                                                    if(self.B[i][l][m]!=[] and self.B[j][n][o]!=[] and self.B[k][p][q]!=[]):
                                                        new_L[i][j][k] += self.B[i][l][m].dot(self.B[j][n][o].dot(self.B[k][p][q].dot(self.O[l][m][n][o][p][q])))   

                                                    #################
                                                        
                                                    for r in range(7):
                                                        if(self.A[i][n]!=[] and self.A[j][o]!=[] and self.A[k][p]!=[] and self.A[l][q]!=[] and self.A[m][r]!=[]):
                                                            new_N[i][j][k][l][m] += self.A[i][n].dot(self.A[j][o].dot(self.A[k][p].dot(self.A[l][q].dot(self.A[m][r].dot(self.N[n][o][p][q][r])))))

                                                        if(self.A[i][m]!=[] and self.A[j][n]!=[] and self.B[k][o][p]!=[] and self.B[l][q][r]!=[]):
                                                            new_M[i][j][k][l] += self.A[i][m].dot(self.A[j][n].dot(self.B[k][o][p].dot(self.B[l][q][r].dot(self.O[m][n][o][p][q][r]))))                                                   
                                                        if(self.A[i][m]!=[] and self.B[j][n][o]!=[] and self.A[k][p]!=[] and self.B[l][q][r]!=[]):
                                                            new_M[i][j][k][l] += self.A[i][m].dot(self.B[j][n][o].dot(self.A[k][p].dot(self.B[l][q][r].dot(self.O[m][n][o][p][q][r]))))                                                   
                                                        if(self.A[i][m]!=[] and self.B[j][n][o]!=[] and self.B[k][p][q]!=[] and self.A[l][r]!=[]):
                                                            new_M[i][j][k][l] += self.A[i][m].dot(self.B[j][n][o].dot(self.B[k][p][q].dot(self.A[l][r].dot(self.O[m][n][o][p][q][r]))))                                                   
                                                        if(self.B[i][m][n]!=[] and self.B[j][o][p]!=[] and self.A[k][q]!=[] and self.A[l][r]!=[]):
                                                            new_M[i][j][k][l] += self.B[i][m][n].dot(self.B[j][o][p].dot(self.A[k][q].dot(self.A[l][r].dot(self.O[m][n][o][p][q][r]))))     

                                                        if(self.B[i][m][n]!=[] and self.A[j][o]!=[] and self.B[k][p][q]!=[] and self.A[l][r]!=[]):
                                                            new_M[i][j][k][l] += self.B[i][m][n].dot(self.A[j][o].dot(self.B[k][p][q].dot(self.A[l][r].dot(self.O[m][n][o][p][q][r]))))                                                   
                                                        if(self.B[i][m][n]!=[] and self.A[j][o]!=[] and self.A[k][p]!=[] and self.B[l][q][r]!=[]):
                                                            new_M[i][j][k][l] += self.B[i][m][n].dot(self.A[j][o].dot(self.A[k][p].dot(self.B[l][q][r].dot(self.O[m][n][o][p][q][r]))))                             

                self.set_Carleman(J_in = new_J, K_in = new_K, L_in= new_L, M_in= new_M, N_in= new_N, O_in= new_O)

        else: # Direct evolution
            new_J = self.J.copy()
            new_J[0] += -self.dt*(self.D[0].dot(self.J[1])+self.D[1].dot(self.J[2])) 
            new_J[1] += -self.dt*(self.D[0].dot(self.J[3])+self.D[1].dot(self.J[4])) 
            new_J[2] += -self.dt*(self.D[0].dot(self.J[5])+self.D[1].dot(self.J[6])) 
            new_J[3] += -self.dt*(self.omega*self.J[3]
                                 -self.omega*self.J[1]*self.J[1]*(2-self.J[0])
                                 -self.cs**2*self.omega*self.J[0]
                                 +self.cs**2*3*self.cs**2*self.D[0].dot(self.J[1])
                                 +self.cs**2*self.D[1].dot(self.J[2]))
            new_J[4] += -self.dt*(self.omega*self.J[4]
                                 -self.omega*self.J[1]*self.J[2]*(2-self.J[0])
                                 +self.cs**2*self.D[1].dot(self.J[1])
                                 +self.cs**2*self.D[0].dot(self.J[2]))
            new_J[5] += -self.dt*(self.omega*self.J[5]
                                 -self.omega*self.J[1]*self.J[2]*(2-self.J[0])
                                 +self.cs**2*self.D[1].dot(self.J[1])
                                 +self.cs**2*self.D[0].dot(self.J[2]))
            new_J[6] += -self.dt*(self.omega*self.J[6]
                                 -self.omega*self.J[2]*self.J[2]*(2-self.J[0])
                                 -self.cs**2*self.omega*self.J[0]
                                 +self.cs**2*self.D[0].dot(self.J[1])
                                 +self.cs**2*3*self.D[1].dot(self.J[2]))
           
    def print_to_file(self, filepath):
        for x in range(self.sites):
            if(x == 0): df = pd.DataFrame([({**{f'x_{d}': int(self.coord(x)[d],self.length) for d in range(2)},**{f'J_{n}':self.J[n][x] for n in range(len(self.J))}})])
            else: df = pd.concat([df,pd.DataFrame([({**{f'x_{d}': int(self.coord(x)[d],self.length) for d in range(2)},**{f'J_{n}':self.J[n][x] for n in range(len(self.J))}})])])
        df.to_csv(filepath, float_format='%.8f',mode= 'a', index=False, header=not os.path.exists(filepath))


    def non_local_evolution(self):

        new_J = []
        new_K = []
        for i in range(7):
            new_J.append(np.zeros(self.sites))
            new_K.append([])
            for j in range(7):
                if(self.A[i][j]!=[]):
                    new_J[i] += self.A[i][j].dot(self.J[j])
                new_K[-1].append([])
                for k in range(7):
                    if(self.B[i][j][k]!=[]):
                        for x in range(self.sites):
                            new_J[i][x] += self.B[i][j][k].toarray()[x][x]*self.K[j][k][x][x]
                    for l in range(7):
                        if(self.A[i][k]!=[] and self.A[j][l]!=[]):
                            print(i,j,k,l)
                            # for x in range(self.sites):
                            #     for y in range(self.sites):
                            #         for x1 in range(self.sites):
                            #             for y1 in range(self.sites):
                            #                 new_K[i][j][x][y] += self.A[i][k].toarray()[x][x1]*self.A[j][l].toarray()[y][y1]*self.K[k][l][x1][y1]           
                            new_K[i][j] = np.einsum('ij,kl,jl->ik',self.A[i][k].toarray(),self.A[j][l].toarray(),self.K[k][l])
        self.set_Carleman(J_in = new_J, K_in = new_K)
                


