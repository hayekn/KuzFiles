import numpy as np
import matplotlib.pyplot as plt 

plt.plot(*np.loadtxt("test.dat",unpack=True, usecols=(0,1)), linewidth=2.0)
plt.show()

#0=TIME
#1=2nd col, etc...