import numpy as np
import matplotlib.pyplot as plt 

'''plt.plot(*np.loadtxt("control.dat",unpack=True, usecols=(0,2)), linewidth=2.0, label='Binge')
plt.plot(*np.loadtxt("control.dat",unpack=True, usecols=(0,3)), linewidth=2.0, label='Stop')
plt.plot(*np.loadtxt("control.dat",unpack=True, usecols=(0,1)), linewidth=2.0, label='Seek')
plt.plot(*np.loadtxt("control.dat",unpack=True, usecols=(0,4)), linewidth=2.0, label='NAc')'''
plt.plot(*np.loadtxt("control.dat",unpack=True, usecols=(0,5)), linewidth=2.0, label='Volume')
plt.legend()
plt.show()

#0=TIME
#1=2nd col, etc...