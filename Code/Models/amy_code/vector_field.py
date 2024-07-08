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
y_traj = [0, 0, 0.1, 0, 0.3, 0.3, 0] 
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

def der_model(t, y):
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



#Plotting the Vectorfield
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
                y = [init[0,z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z], init[6,z]]
                y[n] = x1array[i, j]
                y[m] = x2array[i, j]
                deriv = der_model(t, y)
                dx1dt_array[i, j] = deriv[n]
                dx2dt_array[i, j] = deriv[m]
        ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)
        y = xppaut_model(t, y_traj)['Int']
        traj = np.zeros((7, len(t)))
        for k in np.arange(7):
            traj[k, :] = y[k]  #seek, setp, binge, spike, nac, av, ALCOHOL
        
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
            
            vta = F(Evta*(y0[6] - TOLERANCE*daFACTOR)) 
            #PFC System
            nullcline[0,i] = (F(Eseek * (spTOseek * y0[1] - apsTOseek * y0[5] - seekDRIVE ))) / seekTAU
            nullcline[1,i] = (F(Esetp * (TOLERANCE - y0[6]))) / setpTAU

            #Insular System (sustained+spike)
            nullcline[2,i] = (F(Esustained * (- seekTOsus * y0[0] - sustainedDRIVE ))) / sustainedTAU
            nullcline[3,i] = (F(Espike * (((time[i]-spikeDUR)- seekTOspike * y0[0] -spikeDRIVE )))) / spikeTAU

            insula = (y0[2] + y0[3])/insula_norm
            #Striatal System
            nullcline[4,i]= (F(Enac * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * insula - nacDRIVE ))) / nacTAU
            #Alcohol Tracking
            nullcline[5,i] = (F(Eav * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * insula - avDRIVE ))) / avTAU
              
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


vector_field(y0,y_traj,t,0,2, ['Seek', 'Binge'], 'no')



# def vec_field(y0, y_traj,t, n, m, name, save):
#     #Setting up the Meshgrid and Subplots
#     x1min, x1max, numptsx1, x2min, x2max, numptsx2 = -0.1, 1, 12, -0.1, 1, 12
#     x1list = np.linspace(x1min, x1max, numptsx1)
#     x2list = np.linspace(x2min, x2max, numptsx2)
#     x1array, x2array = np.meshgrid(x1list, x2list)
#     dx1dt_array = np.zeros(x1array.shape)
#     dx2dt_array = np.zeros(x1array.shape)
#     fig = plt.figure(figsize=(12, 10))
#     ax1 = fig.add_subplot(2, 1, 1)     
#     ax2 = fig.add_subplot(2, 2, 3)      
#     ax3 = fig.add_subplot(2, 2, 4)      
#     ax = [ax1, ax2, ax3]
#     init = xppaut_model(t, y0)['Int'] # Solving the system with IC

#     def update(z):
#         ax1.clear()
#         for i in np.arange(numptsx1):
#             for j in np.arange(numptsx2):
#                 y = [init[0, z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z]]
#                 y[n] = x1array[i, j]
#                 y[m] = x2array[i, j]
#                 deriv = derModel(t, y)
#                 dx1dt_array[i, j] = deriv[n]
#                 dx2dt_array[i, j] = deriv[m]
#         ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)
        
#         y = xppaut_model(t, y_traj)['Int']
#         traj = np.zeros((6, len(t)))
#         for k in np.arange(6):
#             traj[k, :] = y[k]  # seek, binge, stop, nac, dls, alc
        
#         setp = 1-F(Esetp*(traj[5,:]-TOLERANCE))

#         ax1.plot(traj[n, 0:z], traj[m, 0:z], color='black', linewidth = '2')
#         ax1.plot(traj[n,z],traj[m,z],marker='o', color = 'red', markersize='10', zorder=10)

#         x1list_fine = np.linspace(x1min, x1max, 250)
#         x2list_fine = np.linspace(x2min, x2max, 250)
#         nullcline = np.zeros((5, 250))
#         for i in np.arange(250):
#             for k in np.arange(6):
#                 y0[k] = traj[k,z]
#             y0[n] = x1list_fine[i]
#             y0[m] = x2list_fine[i]
#             setp_n = np.array(1 - F(Esetp * (y0[5] - TOLERANCE)))
#             vta = F(Evta * (y0[5] - TOLERANCE * daFACTOR))
#             aps = Eaps * (y0[3] + y0[4]) / 2
#             nullcline[0, i] = (F(Eseek * (spTOseek * setp_n - apsTOseek * aps - seekDRIVE))) / seekTAU
#             nullcline[1, i] = (F(Ebinge * (stopTObin * y0[2] - seekTObin * y0[0] - bingeDRIVE))) / bingeTAU
#             nullcline[2, i] = (F(Estop * (binTOstop * y0[1] - spTOstop * setp_n - stopDRIVE))) / stopTAU
#             nullcline[3, i] = (F(Enac * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * y0[1] - nacDRIVE))) / nacTAU
#             nullcline[4, i] = (F(Edls * (-binTOdls * y0[1] - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
        
#         ax1.plot(x1list_fine, nullcline[m, :], 'b-', alpha=0.8, linewidth=1.5)
#         ax1.plot(nullcline[n, :], x2list_fine, 'r-', alpha=0.8, linewidth=1.5)
#         ax1.set_xlabel(name[0], fontweight='bold', fontsize=15, **tfont)
#         ax1.set_ylabel(name[1], fontweight='bold', fontsize=15, **tfont)
#         ax1.set_title('Phase Plane Analysis of '+name[0]+' and '+name[1]+'', **tfont, fontweight = 'bold', fontsize = '15')
        
#         ax2.plot(t[0:z],traj[n,0:z], color = 'red')
#         ax2.plot(t[0:z],traj[m,0:z], color = 'blue')
#         ax2.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
#         ax2.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
#         ax2.set_xlim(0,t[-1])
#         ax2.set_ylim(0,1.05)        
#         ax2.set_title('Insular Activity', **tfont, fontweight = 'bold', fontsize = '15')
#         ax2.legend(name)
    
#         ax3.plot(t[0:z],setp[0:z], color = 'darkgreen', label = 'Setpoint')
#         ax3.set_xlim(0,t[-1])
#         ax3.set_ylim(-0.05,1.05)
#         ax3.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
#         ax3.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
#         ax3.set_title('Setpoint Activity', **tfont, fontweight = 'bold', fontsize = '15')

#         plt.subplots_adjust(hspace=0.3)

#         return ax
#     ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1,repeat=False)

#     if save =='yes':
#         writer = PillowWriter(fps=30)
#         ani.save('/Users/amyrude/Downloads/phase_plane_animation_stable_null.gif', writer=writer)
#     plt.show()
