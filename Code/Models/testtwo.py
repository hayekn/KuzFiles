import numpy as np
from params import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from pythonizedModel import *
from myfunctions import *

y0 = [0, 0.2, 0.2, 0.3, 0, 0]
y_traj = [0, 0.2, 0.2, 0.3, 0, 0]

t = np.linspace(0, 20, 1000)

sub_plots(t, y0)
comb_plots(t,y0)
ind_plots(t,y0)
phase_space(y0,y_traj, t,1,2,['Binge', 'Stop'])
