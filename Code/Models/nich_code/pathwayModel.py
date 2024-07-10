import numpy as np
from pathwayModelParams import *
from datetime import date
import math
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast
from SALib.sample import fast_sampler
import matplotlib.animation as animation

#seek, nseek, seekD1, seekD2, nseekD1, nseekD2, 
#seekGPE, seekGPI, seekSTN, nseekGPE, nseekGPI, nseekSTN,  
#binge, nac, ALCOHOL, aps, setp, vta

y0 = [.2, .2, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      .2, .2, 0, 0, 0, 0]

def F(x):
        return 1 / (1 + np.exp(-x))

def xppaut_model(t, fuzz, nsLEVEL=nsLEVEL, 
                 seekTOseekD1=seekTOseekD1,
                 seekTOseekD2=seekTOseekD2,
                 nseekTOnseekD1=nseekTOnseekD1,
                 nseekTOnseekD2=nseekTOnseekD2):

    def model(t, y):
        seek, nseek, seekD1, seekD2, nseekD1, nseekD2,seekGPE, seekGPI, seekSTN, nseekGPE, nseekGPI, nseekSTN, binge, nac, ALCOHOL, aps, setp, vta = y
        if fuzz:
            noise = np.random.normal(0,.3)
        else:
            noise = 0
            
        ns = nsLEVEL*(F(Ens*(t-nsSTART))+F(Ens*(nsDURATION +nsSTART - t))-1) #check
        cs = 1 - np.heaviside(ALCOHOL - TOLERANCE*daFACTOR, 1)
    
        dvta_dt = -vta + F(Evta*(csTOvta*cs - nsTOvta*ns + vtaDRIVE)) + noise
        dsetp_dt = (-setp + F(Esetp * (ALCOHOL - TOLERANCE)) - F(Esetp * (ALCOHOL - TOLERANCE - setpDELTA))) / setpTAU
        daps_dt = -aps + F(Eaps*nac) + noise
        dseek_dt = (-seek + F(Eseek * (csTOseek * cs - nseekTOseek * nseek - seekGPI + seekDRIVE)) + noise) / seekTAU
        dnseek_dt = (-nseek + F(Enseek * (spTOnseek * setp + nsTOnseek * ns - seekTOnseek * seek - nseekGPI + nseekDRIVE)) + noise) / nseekTAU
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - nsTObin * ns + bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE)) + noise) / nacTAU
        dALCOHOL_dt = nac

        dseekD1_dt = -seekD1 + F(seekTOseekD1 * seek)
        dseekD2_dt = -seekD2 + F(seekTOseekD2 * seek)
        dnseekD1_dt = -nseekD1 + F(nseekTOnseekD1 * nseek)
        dnseekD2_dt = -nseekD2 + F(nseekTOnseekD2 * nseek)

        dseekGPE_dt = -seekGPE + F(-seekD2)
        dseekSTN_dt = -seekSTN + F(-seekGPE)
        dseekGPI_dt = -seekGPI + F(seekSTN - seekD1)

        dnseekGPE_dt = -nseekGPE + F(-nseekD2)
        dnseekSTN_dt = -nseekSTN + F(-nseekGPE)
        dnseekGPI_dt = -nseekGPI + F(nseekSTN - nseekD1)

        return [dseek_dt, dnseek_dt, dseekD1_dt, dseekD2_dt, dnseekD1_dt, dnseekD2_dt, 
                dseekGPE_dt, dseekGPI_dt, dseekSTN_dt, dnseekGPE_dt, dnseekGPI_dt, dnseekSTN_dt,  
                dbinge_dt, dnac_dt, dALCOHOL_dt, daps_dt, dsetp_dt, dvta_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs(time=120, slide=1, quick=1, decay=1, fuzz=False, save=False, show=False):
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    t = np.linspace(0, time, 300)
    y = xppaut_model(t, fuzz, nsLEVEL=0)
    expectedReward = 0

    SD1=seekTOseekD1
    SD2=seekTOseekD2
    ND1=nseekTOnseekD1
    ND2=nseekTOnseekD2
    
    for i in range(slide):
        y = xppaut_model(t, fuzz, nsLEVEL=0, 
                         seekTOseekD1=SD1, seekTOseekD2=SD2,
                         nseekTOnseekD1=ND1, nseekTOnseekD2=ND2)
        rewardError = (y['Int'][14][-1] - expectedReward)/y['Int'][14][-1]
        expectedReward = (expectedReward + y['Int'][14][-1])/2

        if y['Int'][0][-1] > y['Int'][1][-1]:
            SD1=SD1*math.exp(rewardError*quick)
            SD2=SD2*math.exp(-rewardError*quick)
            print("yes")
        else:
            ND1=ND1*math.exp(rewardError*quick)
            ND2=ND2*math.exp(-rewardError*quick)
            print("no")

        SD1 = SD1*math.exp((1-SD1)*decay)
        SD2 = SD2*math.exp((1-SD2)*decay)
        ND1 = ND1*math.exp((1-ND1)*decay)
        ND2 = ND2*math.exp((1-ND2)*decay)
            
    print(ND1, ND2, SD1, SD2)


    cs = 1 - np.heaviside(y['Int'][14] - TOLERANCE*daFACTOR, 1)
    seek = axs[0, 0].plot(t, y['Int'][0], label="Seeking", color = 'midnightblue')[0]
    nseek = axs[0, 0].plot(t, y['Int'][1], label="Not Seeking", color='blue')[0]
    sp = axs[0, 0].plot(t, y['Int'][16], label="Setpoint", color = 'royalblue')[0]
    bin = axs[0, 1].plot(t, y['Int'][12], label="binge", color = 'mediumseagreen')[0]
    nac = axs[0, 2].plot(t, y['Int'][13], label="NAc", color = 'maroon')[0]
    av = axs[0, 2].plot(t, y['Int'][15], '--', label='AV',color='maroon')[0]
    alc = axs[1, 0].plot(t, y['Int'][14], label='Alcohol Vol.', color = 'red')[0]
    da = axs[1, 1].plot(t, y['Int'][17], label="VTA", color = 'lightcoral')[0]
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

    if save:
        plt.savefig("newGraphs"+str(date.today()), dpi=350)
    else:
        plt.show()
    
runGraphs(time=75, slide=100, quick=1, decay=.01)

