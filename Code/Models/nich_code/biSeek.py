import numpy as np
from biSeekParams import *
from datetime import date
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast
from SALib.sample import fast_sampler
import matplotlib.animation as animation

#seek, binge, nac, ALCOHOL, aps, setp, vta, nseek

def F(x):
        return 1 / (1 + np.exp(x))

def xppaut_model(t, fuzz, y0=y0, 
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
vtaDRIVE = vtaDRIVE,
setpTAU = setpTAU,

#SYNAPTIC WEIGHTS
spTOseek=spTOseek,
seekTOnac=seekTOnac,
seekTObin=seekTObin,
binTOnac=binTOnac,
binTObin=binTObin,
vtaTOnac=vtaTOnac,
apsTOseek=apsTOseek,
csTOvta = csTOvta,
csTOseek = csTOseek,

#NEGATIVE STIM
Ens = Ens,
nsLEVEL = nsLEVEL,
nsSTART=nsSTART,
nsDURATION=nsDURATION,
nsTOvta = nsTOvta,
nsTOseek = nsTOseek,
nsTObin = nsTObin,

#EXTRAS
TOLERANCE=TOLERANCE,
daFACTOR=daFACTOR,
decayFac=decayFac
):

    def model(t, y):
        seek, binge, nac, ALCOHOL, aps, setp, vta, nseek = y
        if fuzz:
            noise = np.random.normal(0,.3)
        else:
            noise = 0
            
        ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-nsSTART-nsDURATION))-1)
        cs = 1 - np.heaviside(ALCOHOL - TOLERANCE*daFACTOR, 1)
    
        dvta_dt = -vta + F(Evta*(-csTOvta*cs + nsTOvta*ns - vtaDRIVE)) + noise
        dsetp_dt = (-setp + F(Esetp * (TOLERANCE - ALCOHOL)) - F(Esetp * (TOLERANCE - ALCOHOL + setpDELTA))) / setpTAU
        daps_dt = -aps + F(-Eaps*nac) + noise
        dseek_dt = (-seek + F(Eseek * (-csTOseek * cs + nseekTOseek * nseek - seekDRIVE)) + noise) / seekTAU
        dnseek_dt = (-nseek + F(Enseek * (-spTOnseek * setp - nsTOnseek * ns + seekTOnseek * seek - nseekDRIVE)) + noise) / nseekTAU
        dbinge_dt = (-binge + F(Ebinge * (-seekTObin * seek + nsTObin * ns - bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE)) + noise) / nacTAU
        dALCOHOL_dt = nac

        return [dseek_dt, dbinge_dt, dnac_dt, dALCOHOL_dt, daps_dt, dsetp_dt, dvta_dt, dnseek_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs(time=120, nsSTOP=40, nsDURATION=nsDURATION, frames=100, fuzz=False, nsAnimation=False, save=False):
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    t = np.linspace(0, time, 300)

    if nsAnimation:
        y = xppaut_model(t, fuzz)
        ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-(nsSTART+nsDURATION)))-1)
    else:
        y = xppaut_model(t, fuzz, nsLEVEL=0)
        ns = np.zeros(len(t))

    cs = 1 - np.heaviside(y['Int'][3] - TOLERANCE*daFACTOR, 1)

    seek = axs[0, 0].plot(t, y['Int'][0], label="Seeking", color = 'midnightblue')[0]
    nseek = axs[0, 0].plot(t, y['Int'][7], label="Not Seeking", color='blue')[0]
    sp = axs[0, 0].plot(t, y['Int'][5], label="Setpoint", color = 'royalblue')[0]
    bin = axs[0, 1].plot(t, y['Int'][1], label="binge", color = 'mediumseagreen')[0]
    nac = axs[0, 2].plot(t, y['Int'][2], label="NAc", color = 'maroon')[0]
    av = axs[0, 2].plot(t, y['Int'][4], '--', label='AV',color='maroon')[0]
    alc = axs[1, 0].plot(t, y['Int'][3], label='Alcohol Vol.', color = 'red')[0]
    da = axs[1, 1].plot(t, y['Int'][6], label="VTA", color = 'lightcoral')[0]
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
    axs[1, 0].set_ylim(0, 50)

    def negAnim(frame):
        newNS = nsSTOP*(frame/frames)
        newY = xppaut_model(t, fuzz, nsSTART=newNS)
        
        ns = nsLEVEL*(F(Ens*(newNS-t))+F(Ens*(t-(newNS+nsDURATION)))-1)
        cs = 1 - np.heaviside(newY['Int'][3] - TOLERANCE*daFACTOR, 1)

        seek.set_data(t, newY['Int'][0])
        nseek.set_data(t, newY['Int'][7])
        bin.set_data(t, newY['Int'][1])
        nac.set_data(t, newY['Int'][2])
        av.set_data(t, newY['Int'][4])
        alc.set_data(t, newY['Int'][3])
        sp.set_data(t, newY['Int'][5])
        da.set_data(t, newY['Int'][6])
        neg.set_data(t, ns)
        cond.set_data(t, cs)

        return (seek, nseek, bin, nac, alc, sp, da, neg, cond)
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
    
#runGraphs(75, fuzz=True)
runGraphs(75,fuzz=True,  nsAnimation=True, save=True)

