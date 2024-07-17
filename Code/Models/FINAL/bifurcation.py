from final import *
from scipy.optimize import fsolve

def system(inputs, setpAct):

    return [(-inputs[0] + F(-spTOseek * setpAct + Eseek * (binTOseek * inputs[1] + seekDRIVE))) / seekTAU,
    (-inputs[1] + F(Ebinge * (seekTObin * inputs[0] + bingeDRIVE))) / bingeTAU]

def equiFinder(i, j, setpAct):
    return fsolve(system, [i, j], args=(setpAct))

def seekBif(limit, reso, save=False):
    E = np.linspace(0, limit, reso)

    for i in np.linspace(0, 1, 15): 
        for j in np.linspace(0, 1, 15):
            axs[0].scatter(E, [equiFinder(i, j, setpAct=e)[0] for e in E], c='black', s=.1)

    axs[0].set_xlabel("Setp Activity")
    axs[0].set_ylabel("Equilibrium Value")
    axs[0].set_ylim(0, 1)

def bingeBif(limit, reso, save=False):
    E = np.linspace(0, limit, reso)

    for i in np.linspace(0, 1, 15): 
        for j in np.linspace(0, 1, 15):
            axs[1].scatter(E, [equiFinder(i, j, setpAct=e)[1] for e in E], c='black', s=.1)

    axs[1].set_xlabel("Setp Activity")
    axs[1].set_ylabel("Equilibrium Value")
    axs[1].set_ylim(0, 1)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

seekBif(1, 1000)
bingeBif(1, 1000)
plt.show()