#Assembly
import numpy as np
import matplotlib.pyplot as plt

#Data InFlow
x = np.loadtxt("GridCells.txt")

#Func
plt.pcolor(x)
plt.show()
