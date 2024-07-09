import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import random
# from myfunctions import sub_plots
tfont = {'fontname':'Times New Roman'}

y0 = [0, 0, 0.1, 0, 0.3, 0.3, 0] #seek, setp, binge, spike, nac, av, ALCOHOL
t= np.linspace(0,50,500)  

## Defining the Parameters

#EXCITABILITIES
Esustained = 5
Enac = 1.84
Eav = 1.84
Eaps = 1
Eseek=1
Esetp = 15
Evta = 12
Espike = 5

#TIMESCALESf
seekTAU = 1
sustainedTAU = 1
nacTAU = 1
avTAU =1
spikeTAU=1
setpTAU =1

#DRIVES
seekDRIVE = 0.01
sustainedDRIVE = -0.5
nacDRIVE = -1.6
avDRIVE = -1.4
spikeDRIVE = -4.2

#SYNAPTIC WEIGHTS
spTOseek = 5
seekTOnac = 10
seekTOsus = 3
binTOnac = 1
vtaTOnac = 1
apsTOseek = 1
seekTOspike =0.5

#EXTRAS
TOLERANCE = 20
daFACTOR = 0.1
insula_norm = 1.2
spikeDUR = t[-1]/10

def F(x):
        return 1 / (1 + np.exp(x))
def xppaut_model(t, y0):

    def model(t, y):
        seek, setp, sustained, spike, nac, av, ALCOHOL = y
        #Dopaminergic System
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR)) 
        #PFC System
        dsetp_dt = (-setp + F(Esetp * (TOLERANCE - ALCOHOL))) / setpTAU
        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * av - seekDRIVE ))) / seekTAU
        #Insular System (sustained+spike)
        dsustained_dt = (-sustained + F(Esustained * (- seekTOsus * seek - sustainedDRIVE ))) / sustainedTAU
        dspike_dt = (-spike + F(Espike * (((t-spikeDUR)- seekTOspike * seek -spikeDRIVE )))) / spikeTAU

        insula = (sustained + spike)/insula_norm
        #Striatal System
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * insula - nacDRIVE ))) / nacTAU
        #Alcohol Tracking
        dav_dt = (-av + F(Eav * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * insula - avDRIVE ))) / avTAU
        dALCOHOL_dt = nac
        
        return [dseek_dt, dsetp_dt, dsustained_dt, dspike_dt, dnac_dt, dav_dt, dALCOHOL_dt ]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    insula = (y[2]+y[3])/insula_norm
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'Ins':insula}

def xppaut_model_noise(t, y0):

    def model(t, y):
        seek, setp, sustained, spike, nac, av, ALCOHOL = y
        #Dopaminergic System
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR)) + random.uniform(-0.05,0.05)
        #PFC System
        dsetp_dt = (-setp + F(Esetp * (TOLERANCE - ALCOHOL)) + random.uniform(-0.75,0.75)) / setpTAU
        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * av - seekDRIVE ))+ random.uniform(-0.75,0.75)) / seekTAU
        #Insular System (sustained+spike)
        dsustained_dt = (-sustained + F(Esustained * (- seekTOsus * seek - sustainedDRIVE ))+ random.uniform(-0.75,0.75)) / sustainedTAU
        dspike_dt = (-spike + F(Espike * (((t-spikeDUR)- seekTOspike * seek -spikeDRIVE )))+ random.uniform(-0.75,0.75)) / spikeTAU

        insula = (sustained + spike)/insula_norm
        #Striatal System
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * insula - nacDRIVE ))+ random.uniform(-0.75,0.75)) / nacTAU
        #Alcohol Tracking
        dav_dt = (-av + F(Eav * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * insula - avDRIVE ))+ random.uniform(-0.75,0.75)) / avTAU
        dALCOHOL_dt = nac
        
        return [dseek_dt, dsetp_dt, dsustained_dt, dspike_dt, dnac_dt, dav_dt, dALCOHOL_dt ]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    insula = (y[2]+y[3])/insula_norm
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'Ins':insula}



#Plotting the Data
def sub_plots(t,y0):
    y= xppaut_model(t,y0)
    seek = y['Int'][0]
    setp = y['Int'][1]
    sustained = y['Int'][2]
    spike = y['Int'][3]
    insula = y['Ins']
    nac = y['Int'][4]
    alc = y['Int'][6]
    vta = F(Evta*(alc - TOLERANCE*daFACTOR))+ random.uniform(-0.05,0.05)
    for n in np.arange(len(alc)):
        if alc[n]>=TOLERANCE:
            thresh = t[n] #time at which threshold is reached 
            index = n #index of when threshold is reached
            break
    f, axs = plt.subplots(2, 2, figsize=(10, 8))

    axs[0,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,0].plot(t, seek+setp, label = 'Combined', color = 'lightsteelblue')
    axs[0,0].plot(t, seek, label = 'Seek', color = 'midnightblue')
    axs[0,0].plot(t, setp, label = 'Setpoint', color = 'royalblue')
    axs[0,0].set_title('mPFC Activity', **tfont, fontweight = 'bold', fontsize='14')
    axs[0,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,0].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,0].legend()

    axs[0,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,1].plot(t, sustained/insula_norm, label ='Sustained', color = 'lightgrey')
    axs[0,1].plot(t, spike/insula_norm, label ='Spike', color = 'grey')
    axs[0,1].plot(t, insula, label ='Binge Total', color = 'darkgreen')

    axs[0,1].set_title('Insular Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[0,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,1].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,1].legend()

    axs[1,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,0].plot(t, vta, label = 'DA', color = 'lightcoral')
    axs[1,0].plot(t,nac, label = 'NAc', color = 'maroon')
    axs[1,0].set_title('Subcortical Nuclei Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,0].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[1,0].legend()

    axs[1,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,1].plot(t, alc, color = 'red')
    axs[1,1].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,1].set_ylabel('Volume (mL)',**tfont, fontsize='12')
    axs[1,1].set_xlabel('Time (min)',**tfont, fontsize='12')
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    plt.show()

sub_plots(t,y0)