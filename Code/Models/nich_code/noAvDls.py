import numpy as np
from noAvDlsParams import *
from datetime import date
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast
from SALib.sample import fast_sampler
import matplotlib.animation as animation

#setp, seek, binge, nac, vta, av, ALCOHOL, dls

def F(x):
        return 1 / (1 + np.exp(-x))

def xppaut_model(t, fuzz, nsLEVEL=nsLEVEL, csTOvta=csTOvta):
    y0 = [0, .5, .1, 0.2, 0, 0, 0, EnacMEAN, driveMEAN]
    def model(t, y):
        setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE = y
        if fuzz:
            noise = np.random.normal(0,.3)
        else:
            noise = 0

        ns = nsLEVEL*(F(Ens*(nsSTART-t))+F(Ens*(t-nsSTART-nsDURATION))-1)
        cs = np.heaviside(csDUR-t, 1)

        dEnac_dt = vta/EnacSpeed + EnacDecay*(EnacMEAN-Enac)
        dnacDRIVE_dt =  vta/nacDriveSpeed + NacDriveDecay*(driveMEAN-nacDRIVE)
        
        dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsLEVEL*dls) + setpDRIVE)) + noise) / setpTAU
        dseek_dt = (-seek + F(Eseek * (-spTOseek * setp - nsTOseek * ns + csTOseek * cs + binTOseek * binge + seekDRIVE)) + noise) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - nsTObin * ns + bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE)) + noise) / nacTAU
        ddls_dt = (-dls + F(Edls * (dlsTOdls * dls + csTOdls * cs + dlsDRIVE + noise))) / dlsTAU
        
        dvta_dt = (-vta + F(Evta*(csTOvta * cs - nsTOvta * ns + vtaDRIVE)) + noise) / vtaTAU
        dALCOHOL_dt = (nac+dlsLEVEL*dls)

        return [dsetp_dt, dseek_dt, dbinge_dt, dnac_dt, ddls_dt, dvta_dt, dALCOHOL_dt, dEnac_dt, dnacDRIVE_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    y[4] = dlsLEVEL*y[4]
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs(time=120, fuzz=False, save=False, anim=False):
    fig, axs = plt.subplots(2, 3, figsize=(11, 8))
    t = np.linspace(0, time, 300)
    y = xppaut_model(t, fuzz, nsLEVEL=0)
    
    ns = np.zeros(len(t))
    cs = np.heaviside(csDUR-t, 1)

    sp = axs[0, 0].plot(t, y['Int'][0], label="Setpoint", color = 'royalblue')[0]
    seek = axs[0, 0].plot(t, y['Int'][1], label="Seek", color = 'midnightblue')[0]
    comb = axs[0, 0].plot(t, (y['Int'][0]+y['Int'][1])/2, '--', label="mPFC Average", color = 'lightblue')[0]
    axs[0,0].set_title("mPFC")
    bin = axs[0, 1].plot(t, y['Int'][2], label="Binge", color = 'mediumseagreen')[0]
    axs[0,1].set_title("Insula")
    nac = axs[0, 2].plot(t, y['Int'][3], label="NAc", color = 'maroon')[0]
    dls = axs[0, 2].plot(t, y['Int'][4], label="DLS", color='red')[0]
    sum = axs[0,2].plot(t, y['Int'][3]+y['Int'][4], '--',label="Striatum", color='tomato')[0]
    axs[0, 2].set_title("\"Effort\"")
    da = axs[1, 1].plot(t, y['Int'][5], label="VTA", color = 'lightcoral')[0]
    axs[1,1].set_title("DA")
    alc = axs[1, 0].plot(t, y['Int'][6], label='Alcohol Vol.', color = 'red')[0]
    neg = axs[1, 2].plot(t, ns, label="NS")[0] #ALC OR NEGSTIM OR HEAVISIDE, CHANGE MANUALLY
    cond = axs[1, 2].plot(t, cs, label="CS")[0]
    excNac = axs[1,2].plot(t, y['Int'][7], label="Enac")[0]
    driNac = axs[1, 2].plot(t, np.abs(y['Int'][8]), label='driveNAC')[0]
    axs[1,2].set_title("Conditioned/Negative Stimulus")
    fig.suptitle("Seek <--> Insula; Enac tracks dopamine, decays to baseline (.8); DLS bistable, turned on by CS; no AV variable")
    #Plot formatting
    for i in range(2):
         for j in range(3):
              axs[i, j].legend()
              if i==1 and j==0:
                  continue
              if i==1 and j==2:
                    continue
              axs[i, j].set_ylim(0, 1)
              
    axs[1, 0].set_xlabel('T (min)')
    axs[1, 1].set_xlabel('T (min)')
    axs[1, 2].set_xlabel('T (min)')
    axs[0, 0].set_ylabel('Normalized Activity')
    frames = 100

    def update(frame):
        y = xppaut_model(t, fuzz, nsLEVEL=0, csTOvta=2*(frame/frames))
        sp.set_data(t, y['Int'][0])
        seek.set_data(t, y['Int'][1])
        comb.set_data(t, (y['Int'][0]+y['Int'][1])/2)
        bin.set_data(t, y['Int'][2])
        nac.set_data(t, y['Int'][3])
        dls.set_data(t, y['Int'][4])
        sum.set_data(t, y['Int'][3]+y['Int'][4])
        da.set_data(t, y['Int'][5])
        alc.set_data(t, y['Int'][6])
        excNac.set_data(t, y['Int'][7])
        driNac.set_data(t, np.abs(y['Int'][8]))


        return (sp, seek, nac)
    if anim:
        ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=10)
    
         

    if save and fuzz:
        plt.savefig("newFuzzyGraphs"+str(date.today()), dpi=350)
    elif save and not fuzz:
        plt.savefig("newGraphs"+str(date.today()), dpi=350)
    else:
        plt.show()
    
runGraphs(100, anim=True)

