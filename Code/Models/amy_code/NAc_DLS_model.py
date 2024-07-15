import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import random
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


tfont = {'fontname':'Times New Roman'}
#NOTE: Running this code should produce a subplot with the behavior of all of the neuron groups and a phase-plane analysis animation. 
# The model is defined as binge_model (located below)


#Specify IC and time step
y0 = [0.2, 0, 0.1, 0.1, 0, 0, 0, 0.2] #seek, setp, binge,  nac, dls, ALCOHOL, vta, proxy
y_traj = [0.5, 0, 0.1, 0.2, 0, 0, 0, 0.2] #seek, setp, binge, nac, dls, ALCOHOL, vta, proxy; used in phase plane analysis
t= np.linspace(0,10,10)  

## Defining the Parameters
#EXCITABILITIES
Ebinge = 5
Enac = 1.5
Esetp = 1.5
Eseek= 5
Evta = 3.5
Edls = 1

#TIMESCALES
seekTAU = 1
bingeTAU = 1
nacTAU = 1
setpTAU = 20
vtaTAU = 1
dlsTAU = 1


#DRIVES
seekDRIVE = 1
bingeDRIVE = 1
nacDRIVE = 1.9
setpDRIVE = 1
vtaDRIVE = 1.5
dlsDRIVE = 2.3

#SYNAPTIC WEIGHTS
spTOseek = 8
seekTOnac = 1.5
seekTObin = 1.75
binTOseek = 1.75
binTOnac = 1.5
vtaTOnac = 1.5
csTOseek = 5
csTOvta = 1.5 # modulate this connection to change magnitude of DA peak (1.5 ~ 0.5)
csTOdls = 0.5


#EXTRAS
csDUR = 3
decayFac = 0.001
nacWEIGHT = 0.75

#DA Modulation
proxDECAY = 0.005
EnacSCALE = 2.2
proxTAU = 1.7
nacdrSCALE = 5

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

def binge_model(t, y0, param):
    def model(t, y):
        seek, setp, binge, nac, dls,  ALCOHOL, vta, prox = y

        csTOvta = param
        dprox_dt = vta/proxTAU - (prox * proxDECAY)
        Enac = prox * EnacSCALE
        nacDRIVE = nacdrSCALE/Enac 
        CS = np.heaviside(csDUR-t,0.5) #Conditioned Stimulus
        dseek_dt = (-seek + F(Eseek * (binTOseek * binge + csTOseek * CS - spTOseek * setp - seekDRIVE))) / seekTAU #Seek Activity
        dsetp_dt = (-setp + F(Esetp * (nacWEIGHT * nac + (1-nacWEIGHT) * dls - setpDRIVE - setpDRIVE))) / setpTAU #Alcohol Variable
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - bingeDRIVE))) / bingeTAU #Binge Activity
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge - nacDRIVE))) / nacTAU #NAc Activity
        ddls_dt = (-dls + F(Edls * ( csTOdls * CS - dlsDRIVE)))/ dlsTAU
        dALCOHOL_dt = nac * nacWEIGHT + dls * (1-nacWEIGHT) # Alcohol consumed 
        dvta_dt = (-vta + F(Evta*(csTOvta * CS - vtaDRIVE))) / vtaTAU #VTA activity
       

        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, ddls_dt, dALCOHOL_dt, dvta_dt, dprox_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    CS = np.heaviside(csDUR-t,0.5)
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'CS':CS}

# Some extra functions defined for plotting 
def der_model(t, y , param): #This model is also defined separately, just for convenience. The system of ODEs is also defined in the "binge_model"
        seek, setp, binge, nac, dls,  ALCOHOL, vta, prox = y

        csTOvta = param
        dprox_dt = vta/proxTAU - (prox * proxDECAY)
        Enac = prox * EnacSCALE
        nacDRIVE = nacdrSCALE/Enac 
        CS = np.heaviside(csDUR-t,0.5) #Conditioned Stimulus
        dseek_dt = (-seek + F(Eseek * (binTOseek * binge + csTOseek * CS - spTOseek * setp - seekDRIVE))) / seekTAU #Seek Activity
        dsetp_dt = (-setp + F(Esetp * (nacWEIGHT * nac + (1-nacWEIGHT) * dls - setpDRIVE - setpDRIVE))) / setpTAU #Alcohol Variable
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - bingeDRIVE))) / bingeTAU #Binge Activity
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge - nacDRIVE))) / nacTAU #NAc Activity
        ddls_dt = (-dls + F(Edls * ( csTOdls * CS - dlsDRIVE)))/ dlsTAU
        dALCOHOL_dt = nac * nacWEIGHT + dls * (1-nacWEIGHT) # Alcohol consumed 
        dvta_dt = (-vta + F(Evta*(csTOvta * CS - vtaDRIVE))) / vtaTAU #VTA activity
       

        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, ddls_dt, dALCOHOL_dt, dvta_dt, dprox_dt]

param_array = np.linspace(1,2.5,500)



#Plotting the Data
def sub_plots(t,y0, noise, param):
    if noise == 'yes':  
        y= binge_model_noise(t,y0)
    else:
         y = binge_model(t,y0, param)
    seek = y['Int'][0]
    setp = y['Int'][1]
    binge = y['Int'][2]
    nac = y['Int'][3]
    dls = y['Int'][4]
    alc = y['Int'][5]
    vta = y['Int'][6]
    prox = y['Int'][7]
    CS = y['CS']
    # for n in np.arange(len(alc)):
    #     if alc[n]>=TOLERANCE:
    #         thresh = t[n] #time at which threshold is reached 
    #         index = n #index of when threshold is reached
    #         break
    f, axs = plt.subplots(2, 3, figsize=(12, 8))

    # axs[0,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,0].plot(t, (seek+setp)/2,'--' ,label = 'Combined', color = 'lightsteelblue')
    axs[0,0].plot(t, seek, label = 'Seek', color = 'midnightblue')
    axs[0,0].plot(t, setp, label = 'Setpoint', color = 'royalblue')
    axs[0,0].set_title('mPFC Activity', **tfont, fontweight = 'bold', fontsize='14')
    axs[0,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,0].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,0].legend()

    # axs[0,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,1].plot(t,binge, label ='Binge', color = 'mediumseagreen')
    axs[0,1].set_title('Insular Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[0,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,1].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,1].legend()

    # axs[0,2].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,2].plot(t,nac, label = 'NAc', color = 'maroon')
    axs[0,2].plot(t,dls, label = 'DLS', color = 'red')
    axs[0,2].plot(t, nac * nacWEIGHT + (1-nacWEIGHT) * dls,  '--', label = 'Striatum', color = 'pink')
    axs[0,2].set_title('NAc Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[0,2].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,2].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,2].legend()

    # axs[1,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,0].plot(t, vta, label = 'DA', color = 'firebrick')
    axs[1,0].plot(t, prox, label = 'prox' , color = 'blue')

    axs[1,0].set_title('VTA Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,0].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[1,0].legend()

    # axs[1,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,1].plot(t, CS, color = 'peru', label = 'CS')
    # axs[1,1].plot(t, av, color = 'sienna', label = 'AV')
    axs[1,1].set_title('Conditioned Stimulus and Alc Var',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,1].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[1,1].legend()

    # axs[1,2].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,2].plot(t, alc, color = 'red')
    axs[1,2].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,2].set_ylabel('Volume (mL)',**tfont, fontsize='12')
    axs[1,2].set_xlabel('Time (min)',**tfont, fontsize='12')
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.show()

def sub_plots_ani(t,y0, param_array, save):
   f, axs = plt.subplots(2, 3, figsize=(12, 8))
   def update(z):
        axs[0,0].clear()
        axs[0,1].clear()
        axs[0,2].clear()
        axs[1,0].clear()
        axs[1,1].clear()
        axs[1,2].clear()

        new = param_array[z]
        y = binge_model(t,y0,new)
        seek = y['Int'][0]
        setp = y['Int'][1]
        binge = y['Int'][2]
        nac = y['Int'][3]
        dls = y['Int'][4]
        alc = y['Int'][5]
        vta = y['Int'][6]
        prox = y['Int'][7]
        CS = y['CS'] 
        axs[0,0].plot(t, (seek+setp)/2,'--' ,label = 'Combined', color = 'lightsteelblue')
        axs[0,0].plot(t, seek, label = 'Seek', color = 'midnightblue')
        axs[0,0].plot(t, setp, label = 'Setpoint', color = 'royalblue')
        axs[0,0].set_title('mPFC Activity', **tfont, fontweight = 'bold', fontsize='14')
        axs[0,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[0,0].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[0,0].legend()

        axs[0,1].plot(t,binge, label ='Binge', color = 'mediumseagreen')
        axs[0,1].set_title('Insular Activity',**tfont, fontweight = 'bold', fontsize='14')
        axs[0,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[0,1].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[0,1].legend()

        axs[0,2].plot(t,nac, label = 'NAc', color = 'maroon')
        axs[0,2].plot(t,dls, label = 'DLS', color = 'red')
        axs[0,2].plot(t, nac * nacWEIGHT + (1-nacWEIGHT) * dls,  '--', label = 'Striatum', color = 'pink')
        axs[0,2].set_title('NAc Activity',**tfont, fontweight = 'bold', fontsize='14')
        axs[0,2].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[0,2].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[0,2].legend()
        axs[0,2].set_ylim(0,1)

        axs[1,0].plot(t, vta, label = 'DA', color = 'firebrick')
        axs[1,0].set_ylim(0,1)

        # axs[1,0].plot(t, prox * EnacSCALE, label = 'Enac' , color = 'blue')

        axs[1,0].set_title('VTA Activity',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[1,0].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,0].legend()

        axs[1,1].plot(t, CS, color = 'peru', label = 'CS')
        axs[1,1].set_title('Conditioned Stimulus and Alc Var',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[1,1].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,1].legend()

        axs[1,2].plot(t, alc, color = 'red')
        axs[1,2].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,2].set_ylabel('Volume (mL)',**tfont, fontsize='12')
        axs[1,2].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,2].set_ylim(0,10)

        plt.subplots_adjust(wspace=0.3, hspace=0.3)
        return axs
   ani = animation.FuncAnimation(f, update, frames=len(param_array), interval=1,repeat=False)
   if save =='yes':
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulateNAc.gif', writer=writer)
   plt.show()


def vector_field(y0, y_traj, t, n,m, name, save): 
    #arguments: initial conditions, t, n/m=neuron population numbers (same as IC), name = ['pop one', 'pop two'], save =='yes' if you want to save]
    x1min, x1max, numptsx1, x2min, x2max, numptsx2 = -0.1, 1, 12, -0.1, 1, 12 #setting up the grid space
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
    init = binge_model(t, y0, csTOvta)['Int'] # Solving the system with IC
    
    def update(z):
        #Plotting the Vector Fields
        ax1.clear()
        for i in np.arange(numptsx1):
            for j in np.arange(numptsx2):
                y = [init[0,z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z], init[6,z], init[7,z]]
                y[n] = x1array[i, j]
                y[m] = x2array[i, j]
                deriv = der_model(t, y, csTOvta)
                seek_der = deriv[0]
                dls_der = deriv[4]
                vta_der = deriv[6]
                deriv[0] = seek_der[z]
                deriv[4] = dls_der[z]
                deriv[6] = vta_der[z]
                dx1dt_array[i, j] = deriv[n]
                dx2dt_array[i, j] = deriv[m]
        ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)
        y = binge_model(t, y_traj, csTOvta)['Int']
        traj = np.zeros((7, len(t)))
        for k in np.arange(7):
            traj[k, :] = y[k]  #seek, setp, binge,  nac, dls, ALCOHOL, vta
        
        ax1.plot(traj[n, 0:z], traj[m, 0:z], color='black', linewidth = '2')
        ax1.plot(traj[n,z],traj[m,z],marker='o', color = 'red', markersize='10', zorder=10)

        x1list_fine = np.linspace(x1min, x1max, 250)
        x2list_fine = np.linspace(x2min, x2max, 250)
        nullcline = np.zeros((7, 250)) #seek, setp, binge,  nac, dls, ALCOHOL, vta
        for i in np.arange(250):
            for k in np.arange(7):
                y0[k] = traj[k,z]
            y0[n] = x1list_fine[i]
            y0[m] = x2list_fine[i]
            #Solving for the nullclines at each time step
            CS = np.heaviside(csDUR-t,0.5)[z] #Conditioned Stimulus
            nullcline[0,i] = (F(Eseek * (binTOseek * y0[2] + csTOseek * CS - spTOseek * y0[1] - seekDRIVE))) / seekTAU  #Seek Activity
            nullcline[1,i] =  (F(Esetp * (nacWEIGHT * y0[3] + (1-nacWEIGHT) * y0[4] - setpDRIVE))) / setpTAU #Setpoint Activity       
            nullcline[2,i] = (F(Ebinge * (seekTObin * y0[0] - bingeDRIVE))) / bingeTAU #Binge Activity
            nullcline[3,i] = (F(Enac * (vtaTOnac * y0[6] + seekTOnac * y0[0] + binTOnac * y0[2] - nacDRIVE))) / nacTAU  #NAc Activity
            nullcline[4,i] = (F(Edls * (csTOdls * CS - dlsDRIVE)))/ dlsTAU #DLS Activity
            # nullcline[5,i] =  # Alcohol consumed 
            nullcline[6,i] = F(Evta*(csTOvta * CS - vtaDRIVE)) #VTA activity    
            
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
        ax3.set_ylim(-0.05,0.3)
        ax3.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
        ax3.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
        ax3.set_title('Setpoint Activity', **tfont, fontweight = 'bold', fontsize = '15')
        return ax
    
    ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1,repeat=False)
    if save =='yes':
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/seek_binge_phaseplane.gif', writer=writer)
    plt.show()

# sub_plots(t, y0, 'no', csTOvta)


# for i in np.arange(numptsx1):
#             for j in np.arange(numptsx2):
#                 y = [init[0,z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z], init[6,z], init[7,z]]
#                 y[n] = x1array[i, j]
#                 y[m] = x2array[i, j]
#                 deriv = der_model(t, y, csTOvta)
#                 seek_der = deriv[0]
#                 dls_der = deriv[4]
#                 vta_der = deriv[6]
#                 deriv[0] = seek_der[z]
#                 deriv[4] = dls_der[z]
#                 deriv[6] = vta_der[z]
#                 dx1dt_array[i, j] = deriv[n]
#                 dx2dt_array[i, j] = deriv[m]
#         ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)


# vector_field(y0,y_traj, t, 0, 2, ['Seek', 'Binge'], 'no')
# sub_plots_ani(t,y0,param_array, 'no')



time = np.linspace(0,50,3)
def td_vect(t,y0, time):

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    init = binge_model(t, y0, csTOvta)['Int'] # Solving the system with IC
    for k in np.arange(len(time)):
        seek_array, binge_array , z = np.meshgrid(np.arange(0, 1.2, 0.2), np.arange(0, 1.2, 0.2), time[k])
        dseek_array = np.zeros(seek_array.shape)
        dbinge_array = np.zeros(binge_array.shape)
        for i in np.arange(len(seek_array)):
            for j in np.arange(len(binge_array)):
                y = [init[0, k], init[1, k], init[2, k], init[3, k], init[4,k], init[5, k], init[6, k], init[7, k]]
                y[0] = seek_array[i, j,0]
                y[2] = binge_array[i, j,0]
                deriv = der_model(t, y, csTOvta)
                seek_d = deriv[0]
                seek_der = seek_d[k]
                dseek_array[i,j] = seek_d[k]
                dbinge_array[i,j] = deriv[2]
        ax.quiver(seek_array, binge_array, time[k], dseek_array, dbinge_array,0,length = 0.3)
    plt.show()
    return ax

td_vect(t,y0, time)  
   
    # seek_array, binge_array, time = np.meshgrid(np.arange(0, 1, 0.2),
    #                   np.arange(0, 1, 0.2),
    #                   np.arange(0, 50, 5))
    # deriv = binge_model(time, y0, csTOvta)['Der'] # Solving the system with IC
    # seek_der = np.zeros(len(seek_array), len(time))
    # binge_der = np.zeros(len(binge_array), len(time))

    # for k in len(time):
    #     deriv_t = deriv[k]
    #     seek_der[:,k] = deriv_t[0]
    #     binge_der[:,k] = deriv_t[2]
    #     t_der = 0

    



# #     for k in len(time):  
# #         for i in np.arange(seek_array):
# #             y = [init[0, time[k]], init[1, time[k]], init[2, time[k]], init[3, time[k]], init[4,time[k]], init[5, time[k]], init[6,time[k]], init[7,time[k]]]
# #             y[0] = seek_array[i]
# #             y[2] = binge_array[i]
# #             deriv = der_model(t, y, csTOvta)
# #             seek_der = deriv[0]
# #             binge_der = deriv[2]
# #             dls_der = deriv[4]
# #             vta_der = deriv[6]
# #             deriv[0] = seek_der[z]
# #             deriv[4] = dls_der[z]
# #             deriv[6] = vta_der[z]

# #             dx1dt_array[i, j] = deriv[n]
# #             dx2dt_array[i, j] = deriv[m]
         
    
#     dseek = 
#     dbinge =
#     dtime = 0
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')     

#     init = binge_model(t, y0, csTOvta)['Int'] # Solving the system with IC
    



    # def update(z):
    #     #Plotting the Vector Fields
    #     ax1.clear()
    #     for i in np.arange(numptsx1):
    #         for j in np.arange(numptsx2):
    #             y = [init[0,z], init[1, z], init[2, z], init[3, z], init[4, z], init[5, z], init[6,z], init[7,z]]
    #             y[n] = x1array[i, j]
    #             y[m] = x2array[i, j]
    #             deriv = der_model(t, y, csTOvta)
    #             seek_der = deriv[0]
    #             dls_der = deriv[4]
    #             vta_der = deriv[6]
    #             deriv[0] = seek_der[z]
    #             deriv[4] = dls_der[z]
    #             deriv[6] = vta_der[z]
    #             dx1dt_array[i, j] = deriv[n]
    #             dx2dt_array[i, j] = deriv[m]
    #     ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.003)
    #     y = binge_model(t, y_traj, csTOvta)['Int']
    #     traj = np.zeros((7, len(t)))
    #     for k in np.arange(7):
    #         traj[k, :] = y[k]  #seek, setp, binge,  nac, dls, ALCOHOL, vta
        
    #     ax1.plot(traj[n, 0:z], traj[m, 0:z], color='black', linewidth = '2')
    #     ax1.plot(traj[n,z],traj[m,z],marker='o', color = 'red', markersize='10', zorder=10)

    #     x1list_fine = np.linspace(x1min, x1max, 250)
    #     x2list_fine = np.linspace(x2min, x2max, 250)
    #     nullcline = np.zeros((7, 250)) #seek, setp, binge,  nac, dls, ALCOHOL, vta
    #     for i in np.arange(250):
    #         for k in np.arange(7):
    #             y0[k] = traj[k,z]
    #         y0[n] = x1list_fine[i]
    #         y0[m] = x2list_fine[i]
    #         #Solving for the nullclines at each time step
    #         CS = np.heaviside(csDUR-t,0.5)[z] #Conditioned Stimulus
    #         nullcline[0,i] = (F(Eseek * (binTOseek * y0[2] + csTOseek * CS - spTOseek * y0[1] - seekDRIVE))) / seekTAU  #Seek Activity
    #         nullcline[1,i] =  (F(Esetp * (nacWEIGHT * y0[3] + (1-nacWEIGHT) * y0[4] - setpDRIVE))) / setpTAU #Setpoint Activity       
    #         nullcline[2,i] = (F(Ebinge * (seekTObin * y0[0] - bingeDRIVE))) / bingeTAU #Binge Activity
    #         nullcline[3,i] = (F(Enac * (vtaTOnac * y0[6] + seekTOnac * y0[0] + binTOnac * y0[2] - nacDRIVE))) / nacTAU  #NAc Activity
    #         nullcline[4,i] = (F(Edls * (csTOdls * CS - dlsDRIVE)))/ dlsTAU #DLS Activity
    #         # nullcline[5,i] =  # Alcohol consumed 
    #         nullcline[6,i] = F(Evta*(csTOvta * CS - vtaDRIVE)) #VTA activity    
            
    #     ax1.plot(x1list_fine, nullcline[m, :], 'b-', alpha=0.8, linewidth=1.5)
    #     ax1.plot(nullcline[n, :], x2list_fine, 'r-', alpha=0.8, linewidth=1.5)
    #     ax1.set_xlabel(name[0], fontweight='bold', fontsize=15, **tfont)
    #     ax1.set_ylabel(name[1], fontweight='bold', fontsize=15, **tfont)
    #     ax1.set_title('Phase Plane Analysis of '+name[0]+' and '+name[1]+'', **tfont, fontweight = 'bold', fontsize = '15')

    #     ax2.plot(t[0:z],traj[n,0:z], color = 'red')
    #     ax2.plot(t[0:z],traj[m,0:z], color = 'blue')
    #     ax2.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
    #     ax2.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
    #     ax2.set_xlim(0,t[-1])
    #     ax2.set_ylim(0,1.05)        
    #     ax2.set_title('Neuron Activity', **tfont, fontweight = 'bold', fontsize = '15')
    #     ax2.legend(name)
    
    #     ax3.plot(t[0:z],traj[1,0:z], color = 'darkgreen', label = 'Setpoint')
    #     ax3.set_xlim(0,t[-1])
    #     ax3.set_ylim(-0.05,0.3)
    #     ax3.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
    #     ax3.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
    #     ax3.set_title('Setpoint Activity', **tfont, fontweight = 'bold', fontsize = '15')
    #     return ax
    
    # ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1,repeat=False)
    # if save =='yes':
    #     writer = PillowWriter(fps=30)
    #     ani.save('/Users/amyrude/Downloads/seek_binge_phaseplane.gif', writer=writer)
    # plt.show()