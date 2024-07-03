import numpy as np
from params import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from pythonizedModel import *
from myfunctions import *
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

tfont = {'fontname':'Times New Roman'}

y0 = [0, 0.2, 0.2, 0.3, 0, 0]
y_traj = [0, 0.2, 0.2, 0.3, 0, 0]
t= np.linspace(0,50,1000)  



vec_field(y0, y_traj, t, 2, 1, ['Stop', 'Binge'], 'no')
