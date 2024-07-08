import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import random
# from myfunctions import sub_plots
tfont = {'fontname':'Times New Roman'}

y0 = [0, 0, 0.1, 0.3, 0.3, 0] #seek, setp, binge,  nac, av, ALCOHOL
y_traj = [0, 0, 0.1, 0.3, 0.3, 0] #seek, setp, binge, nac, av, ALCOHOL

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
def der_model(t, y):
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

def vector_field(y0,y_traj, t, n,m, name, save):
    x1min, x1max, numptsx1, x2min, x2max, numptsx2 = -0.1, 1, 12, -0.1, 1, 12
    x1list = np.linspace(x1min, x1max, numptsx1)
    x2list = np.linspace(x2min, x2max, numptsx2)
    x1array, x2array = np.meshgrid(x1list, x2list)
    dx1dt_array = np.zeros(x1array.shape)
    dx2dt_array = np.zeros(x1array.shape)
    fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(2, 1, 1)     
    ax2 = fig.add_subplot(2, 2, 3)      
    ax3 = fig.add_subplot(2, 2, 4)      
    ax = [ax1, ax2, ax3]
    init = xppaut_model(t, y0)['Int'] # Solving the system with IC
    
    def update(z):
        #Plotting the Vector Fields
        ax1.clear()
        for i in np.arange(numptsx1):
            for j in np.arange(numptsx2):
                y = [init[0,z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z]]
                y[n] = x1array[i, j]
                y[m] = x2array[i, j]
                deriv = der_model(t, y)
                dx1dt_array[i, j] = deriv[n]
                dx2dt_array[i, j] = deriv[m]
        ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)
        y = xppaut_model(t, y_traj)['Int']
        traj = np.zeros((6, len(t)))
        for k in np.arange(6):
            traj[k, :] = y[k]  #seek, setp, binge, nac, av, ALCOHOL
        
        ax1.plot(traj[n, 0:z], traj[m, 0:z], color='black', linewidth = '2')
        ax1.plot(traj[n,z],traj[m,z],marker='o', color = 'red', markersize='10', zorder=10)

        x1list_fine = np.linspace(x1min, x1max, 250)
        x2list_fine = np.linspace(x2min, x2max, 250)
        time = np.linspace(0,50,250)
        nullcline = np.zeros((6, 250)) # seek, setp, sustained, spike, nac, av
        for i in np.arange(250):
            for k in np.arange(6):
                y0[k] = traj[k,z]
            y0[n] = x1list_fine[i]
            y0[m] = x2list_fine[i]
            
            vta = F(Evta*(y0[5] - TOLERANCE*daFACTOR)) 
            #PFC System
            nullcline[0,i] = (F(Eseek * (spTOseek * y0[1] - apsTOseek * y0[4] - seekDRIVE ))) / seekTAU
            nullcline[1,i] = (F(Esetp * (TOLERANCE - y0[5]))) / setpTAU

            #Insular System (sustained+spike)
            nullcline[2,i] = ( F(Ebinge * (-vtaTObin*vta - seekTObin * y0[0] - bingeDRIVE))) / bingeTAU

            #Striatal System
            nullcline[3,i]= (F(Enac * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * y0[2] - nacDRIVE ))) / nacTAU
            #Alcohol Tracking
            nullcline[4,i] = (F(Eav * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * y0[2] - avDRIVE ))) / avTAU
              
        ax1.plot(x1list_fine, nullcline[m, :], 'b-', alpha=0.8, linewidth=1.5)
        ax1.plot(nullcline[n, :], x2list_fine, 'r-', alpha=0.8, linewidth=1.5)
        ax1.set_xlabel(name[0], fontweight='bold', fontsize=15, **tfont)
        ax1.set_ylabel(name[1], fontweight='bold', fontsize=15, **tfont)
        ax1.set_title('Phase Plane Analysis of '+name[0]+' and '+name[1]+'', **tfont, fontweight = 'bold', fontsize = '15')

        ax2.plot(t[0:z],traj[n,0:z], color = 'red')
        ax2.plot(t[0:z],traj[m,0:z], color = 'blue')
        ax2.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
        ax2.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
        ax2.set_xlim(0,t[-1])
        ax2.set_ylim(0,1.05)        
        ax2.set_title('Neuron Activity', **tfont, fontweight = 'bold', fontsize = '15')
        ax2.legend(name)
    
        ax3.plot(t[0:z],traj[1,0:z], color = 'darkgreen', label = 'Setpoint')
        ax3.set_xlim(0,t[-1])
        ax3.set_ylim(-0.05,1.05)
        ax3.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
        ax3.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
        ax3.set_title('Setpoint Activity', **tfont, fontweight = 'bold', fontsize = '15')
        return ax
    
    ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1,repeat=False)
    if save =='yes':
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/seek_binge_phaseplane.gif', writer=writer)
    plt.show()


vector_field(y0,y_traj,t,0,2, ['Seek', 'Binge'], 'yes')