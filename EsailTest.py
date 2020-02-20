#Assembly
import numpy as np
import os

#Var
K = 1.38064852e-23
Mp = 1.6726219e-27

#Data Inflow
x,y = np.loadtxt("vtFile.txt",unpack = True)

#Func
sigmax = np.var(x)
sigmay = np.var(y)

Tx = (Mp*sigmax/K)*131.293
print(Tx)