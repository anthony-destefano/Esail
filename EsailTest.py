#Assembly
import numpy as np

#Var
vxtdata = []
vytdata = []

#Data Inflow
f = open("vt_test.txt", "r")
line = f.readline()
cnt = 1
while line:
    x,y = line.strip().split(" ")
    vxtdata.append(x)
    vytdata.append(y)
    line = f.readline()
    cnt += 1

#Func
sigmax = np.var(vxtdata)
sigmay = np.var(vytdata)
