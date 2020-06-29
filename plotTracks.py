#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
# https://matplotlib.org/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
from matplotlib import ticker, cm

filename = sys.argv[1]

n, x, y = np.loadtxt(filename, unpack=True)

plt.scatter(x,y, s=1)
plt.show()