import numpy as np
from params import *
from datetime import date
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast
from SALib.sample import fast_sampler
import matplotlib.animation as animation

#seek, binge, nac, ALCOHOL, aps, setp

def F(x):
        return 1 / (1 + np.exp(x))

def xppaut_model(t, y0=y0, 
Ebinge=Ebinge,
Enac=Enac,
Eaps=Eaps,
Eseek=Eseek,
Esetp=Esetp,
Evta=Evta,

#TIMESCALES
seekTAU=seekTAU,
bingeTAU=bingeTAU,
nacTAU=nacTAU,

#DRIVES
seekDRIVE=seekDRIVE,
bingeDRIVE=bingeDRIVE,
nacDRIVE=nacDRIVE,

#SYNAPTIC WEIGHTS
spTOseek=spTOseek,
seekTOnac=seekTOnac,
seekTObin=seekTObin,
binTOnac=binTOnac,
binTObin=binTObin,
vtaTOnac=vtaTOnac,
apsTOseek=apsTOseek,

#NEGATIVE STIM
Ens = Ens,
nsLEVEL = nsLEVEL,
nsSTART=nsSTART,
nsDURATION=nsDURATION,
nsTOvta = nsTOvta,

#EXTRAS
TOLERANCE=TOLERANCE,
daFACTOR=daFACTOR
):

    def model(t, y):
        seek, binge, nac, ALCOHOL, aps, setp= y

        #setp = F(Esetp*((ALCOHOL - TOLERANCE)))+ F(Esetp*(TOLERANCE - ALCOHOL - spDURATION))-1 #SPECIFIED DURATION SETPOINT
        noise = np.random.normal(0, .3)
        noise = 0

        ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-nsSTART-nsDURATION))-1)
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR + nsTOvta*ns))
        
        dsetp_dt = -setp + np.exp(-decayFac*t)*F(Esetp * (TOLERANCE - ALCOHOL)) + noise
        daps_dt = -aps + F(-Eaps*nac) + noise
        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps + nsTOseek * ns - seekTOseek * seek - seekDRIVE)) + noise) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (-binTObin*binge - seekTObin * seek + spTObin * setp + nsTObinge * ns - vtaTObinge*vta - bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE)) + noise) / nacTAU
        dALCOHOL_dt = nac

        '''for d in [dsetp_dt, daps_dt, dseek_dt, dbinge_dt, dnac_dt]:
            d = d + gaussian_noise'''

        return [dseek_dt, dbinge_dt, dnac_dt, dALCOHOL_dt, daps_dt, dsetp_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}


def runGraphs(time=120, nsSTOP=50, nsDURATION=nsDURATION, frames=100, nsAnimation=False, save=False):
    global fig, axs
    fig, axs = plt.subplots(2, 3, figsize=(10, 6))
    global y0
    y0 = [0, .2, 0.3, 0, 0, 0]
    global t
    t = np.linspace(0, time, 300)
    if nsAnimation:
        y = xppaut_model(t, y0=y0, nsLEVEL=nsLEVEL)
    else:
        y = xppaut_model(t, y0=y0, nsLEVEL=0)
    ALCOHOL = y['Int'][3]
    #setp = np.exp(-decayFac*t)(-setp - F(Esetp * (TOLERANCE - ALCOHOL)))
    ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-(nsSTART+nsDURATION)))-1)
    vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR + ns))
    
    #alcAx = fig.add_subplot(3, 1, 3)
    #alc = alcAx.plot(t, ALCOHOL, label='Consumption')[0]

    global seek
    seek = axs[0, 0].plot(t, y['Int'][0], label="seek")[0]
    axs[0, 0].set_ylabel('Normalized Activity')
    global bin
    bin = axs[0, 1].plot(t, y['Int'][1], label="binge")[0]
    global nac
    nac = axs[0, 2].plot(t, y['Int'][2], label="nac")[0]
    global sp
    sp = axs[1, 0].plot(t, y['Int'][-1], label="setp")[0]
    axs[1, 0].set_ylabel('Normalized Activity')
    global da
    da = axs[1, 1].plot(t, vta, label="vta")[0]
    global neg
    neg = axs[1, 2].plot(t, ALCOHOL, label="ALCOHOL")[0] #ALC OR NEGSTIM, CHANGE MANUALLY
    
    for i in range(2):
         for j in range(3):
              if i==1 and j==2:
                  continue
              axs[i, j].set_xlabel('T')
              axs[i, j].legend()
              #axs[i, j].set_ylim(0, 1)

    axs[1, 2].set_xlabel('T')
    axs[1, 2].legend()

    def negAnim(frame):
        newNS = nsSTOP*(frame/frames)
        newY = xppaut_model(t, y0=y0, nsSTART=newNS)
        ALCOHOL = newY['Int'][3]
        
        ns = nsLEVEL*(F(Ens*(newNS-t))+F(Ens*(t-(newNS+nsDURATION)))-1)
        #setp = (1/(decayFac * t+1))*F(Esetp * (TOLERANCE - ALCOHOL))
        vta = F(Evta*(newY['Int'][3] - TOLERANCE*daFACTOR + ns))
        
        seek.set_data(t, newY['Int'][0])
        bin.set_data(t, newY['Int'][1])
        nac.set_data(t, newY['Int'][2])
        sp.set_data(t, newY['Int'][-1])
        da.set_data(t, vta)
        neg.set_data(t, ALCOHOL) #ALC OR NEGSTIM, CHANGE MANUALLY

        return (seek, bin, nac, sp, da, neg)
    if nsAnimation:
        ani = animation.FuncAnimation(fig=fig, func=negAnim, frames=frames, interval=10)
        if save:
            FFwriter = animation.FFMpegWriter(fps=30)
            ani.save('negStim'+str(date.today())+'.mp4', writer = FFwriter)
            plt.close()
        else:
            plt.show()
    elif save:
        plt.savefig("newGraphs"+str(date.today()), dpi=350)
    plt.show()
    
runGraphs(100)
#runGraphs(nsAnimation=True)

