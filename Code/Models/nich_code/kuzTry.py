import numpy as np
from kuzTryParams import *
from datetime import date
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast
from SALib.sample import fast_sampler
import matplotlib.animation as animation

#setp, seek, binge, nac, vta, av, ALCOHOL

def F(x):
        return 1 / (1 + np.exp(-x))

def xppaut_model(t, fuzz, nsLEVEL=nsLEVEL):

    def model(t, y):
        setp, seek, binge, nac, vta, av, ALCOHOL = y
        if fuzz:
            noise = np.random.normal(0,.3)
        else:
            noise = 0

        ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-nsSTART-nsDURATION))-1)
        cs = np.heaviside(csDUR-t, 1)

        dsetp_dt = (-setp + np.exp(-decayFac*t)*F(Esetp * (avTOsetp * av + setpDRIVE)) + noise) / setpTAU
        dseek_dt = (-seek + F(Eseek * (-spTOseek * setp + avTOseek * av - nsTOseek * ns + csTOseek * cs + seekTOseek * seek + seekDRIVE)) + noise) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - nsTObin * ns + binTObin * binge + bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE)) + noise) / nacTAU
        
        dvta_dt = (-vta + F(Evta*(csTOvta * cs - nsTOvta * ns + vtaDRIVE)) + noise) / vtaTAU
        dav_dt = (-av + F(Eav * (nacTOav*nac + avDRIVE)) + noise) / avTAU
        dALCOHOL_dt = nac

        return [dsetp_dt, dseek_dt, dbinge_dt, dnac_dt, dvta_dt, dav_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs(time=120, fuzz=False, save=False):
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    t = np.linspace(0, time, 300)
    y = xppaut_model(t, fuzz, nsLEVEL=0)
    
    ns = np.zeros(len(t))
    cs = np.heaviside(csDUR-t, 1)

    sp = axs[0, 0].plot(t, y['Int'][0], label="Setpoint", color = 'royalblue')[0]
    seek = axs[0, 0].plot(t, y['Int'][1], label="Seek", color = 'midnightblue')[0]
    bin = axs[0, 1].plot(t, y['Int'][2], label="binge", color = 'mediumseagreen')[0]
    nac = axs[0, 2].plot(t, y['Int'][3], label="NAc", color = 'maroon')[0]
    da = axs[1, 1].plot(t, y['Int'][4], label="VTA", color = 'lightcoral')[0]
    av = axs[0, 2].plot(t, y['Int'][5], '--', label='AV',color='maroon')[0]
    alc = axs[1, 0].plot(t, y['Int'][6], label='Alcohol Vol.', color = 'red')[0]
    neg = axs[1, 2].plot(t, ns, label="NS")[0] #ALC OR NEGSTIM OR HEAVISIDE, CHANGE MANUALLY
    cond = axs[1, 2].plot(t, cs, label="CS")[0]

    #Plot formatting
    for i in range(2):
         for j in range(3):
              if i==1 and j==0:
                  continue
              axs[i, j].set_xlabel('T')
              axs[i, j].legend()
              axs[i, j].set_ylim(0, 1)
    axs[1, 0].set_xlabel('T')
    axs[1, 0].legend()
    axs[0, 0].set_ylabel('Normalized Activity')
    axs[1, 0].set_ylabel('Normalized Activity')

    if save:
        plt.savefig("newGraphs"+str(date.today()), dpi=350)
    else:
        plt.show()
    
runGraphs(75)

