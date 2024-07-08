import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import random
# from myfunctions import sub_plots
tfont = {'fontname':'Times New Roman'}

y0 = [0, 0, 0.1, 0.3, 0.3, 0] #seek, setp, binge, spike, nac, av, ALCOHOL
t= np.linspace(0,50,500)  

## Defining the Parameters

#EXCITABILITIES
Ebinge = 8
Enac = 1.84
Eav = 1.84
Eaps = 1
Eseek=1
Esetp = 15
Evta = 12

#TIMESCALES
seekTAU = 1
bingeTAU = 1
nacTAU = 1
avTAU =1
setpTAU = 1

#DRIVES
seekDRIVE = 0.01
bingeDRIVE = -1.25
nacDRIVE = -1.4
avDRIVE = -1.4

#SYNAPTIC WEIGHTS
spTOseek = 5
seekTOnac = 10
seekTObin = 2.1
binTOnac = 1
vtaTOnac = 1
apsTOseek = 1
vtaTObin = 10

#EXTRAS
TOLERANCE = 20
daFACTOR = 0.1

def F(x):
        return 1 / (1 + np.exp(x))

def xppaut_model(t, y0):

    def model(t, y):
        seek, setp, binge, nac, av, ALCOHOL = y
        #Dopaminergic System
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
        #PFC System
        dsetp_dt = (-setp + F(Esetp * (TOLERANCE - ALCOHOL))) / setpTAU
        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * av - seekDRIVE))) / seekTAU
        #Insular System
        dbinge_dt = (-binge + F(Ebinge * (-vtaTObin*vta - seekTObin * seek - bingeDRIVE))) / bingeTAU
        #Striatal System
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
        #Alcohol Tracking
        dav_dt = (-av + F(Eav * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - avDRIVE))) / avTAU
        dALCOHOL_dt = nac
        
        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, dav_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def xppaut_model_noise(t, y0):

    def model(t, y):
        seek, setp, binge, nac, av, ALCOHOL = y
        #Dopaminergic System
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
        #PFC System
        dsetp_dt = (-setp + F(Esetp * (TOLERANCE - ALCOHOL))+ random.uniform(-0.75,0.75)) / setpTAU
        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * av - seekDRIVE))+ random.uniform(-0.75,0.75)) / seekTAU
        #Insular System
        dbinge_dt = (-binge + F(Ebinge * (-vtaTObin*vta - seekTObin * seek - bingeDRIVE))+ random.uniform(-0.25,0.25)) / bingeTAU
        #Striatal System
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))+ random.uniform(-0.75,0.75)) / nacTAU
        #Alcohol Tracking
        dav_dt = (-av + F(Eav * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - avDRIVE))) / avTAU
        dALCOHOL_dt = nac
        
        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, dav_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

#Plotting the Data
def sub_plots(t,y0):
    y= xppaut_model_noise(t,y0)
    seek = y['Int'][0]
    setp = y['Int'][1]
    binge = y['Int'][2]
    nac = y['Int'][3]
    alc = y['Int'][5]
    vta = F(Evta*(alc - TOLERANCE*daFACTOR))+ random.uniform(-0.08,0.08)
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
    axs[0,1].plot(t,binge, label ='Binge', color = 'mediumseagreen')
   

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