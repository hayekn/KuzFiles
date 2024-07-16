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
y0 = [0.2, 0, 0.1, 0.1, 0, 0, 0, 0.6, -7] #seek, setp, binge,  nac, dls, ALCOHOL, vta, proxy
y_traj = [0.5, 0, 0.1, 0.2, 0, 0, 0, 0.2] #seek, setp, binge, nac, dls, ALCOHOL, vta, proxy; used in phase plane analysis
t= np.linspace(0,30,500)  

## Defining the Parameters
#EXCITABILITIES
Ebinge = 5
Esetp = 1.5
Eseek= 5
Evta = 2
Edls = 3

#TIMESCALES
seekTAU = 1
bingeTAU = 1
nacTAU = 1
setpTAU = 30
vtaTAU = 1
dlsTAU = 1


#DRIVES
seekDRIVE = 1
bingeDRIVE = 1
setpDRIVE = 1
vtaDRIVE = 3.5
dlsDRIVE = 1.4

#SYNAPTIC WEIGHTS
spTOseek = 8
seekTOnac = 2.5
seekTObin = 1.75
binTOseek = 1.75
binTOnac = 2.5
vtaTOnac = 2.5
csTOseek = 5
csTOvta = 4.2 # modulate this connection to change magnitude of DA peak (3 to 5.5)
csTOdls = 2


#EXTRAS
csDUR = 3
decayFac = 0.001
nacWEIGHT = 0.75

#DA Modulation
EnacDECAY = 0.003
nacdriveDECAY = 0.003
EnacTAU = 0.8
nacdriveTAU = 0.8
EnacMEAN = 0.6
driveMEAN = -7
nacdrSCALE = 5
dlsSCALE = 0.1

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

def binge_model(t, y0, param):
    def model(t, y):
        seek, setp, binge, nac, dls,  ALCOHOL, vta, Enac, nacDRIVE = y

        csTOvta = param

        dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN - Enac)
        dnacDRIVE_dt = vta/nacdriveTAU + nacdriveDECAY*(driveMEAN - nacDRIVE)

        CS = np.heaviside(csDUR-t, 0.5) #Conditioned Stimulus
        dseek_dt = (-seek + F(Eseek * (binTOseek * binge + csTOseek * CS - spTOseek * setp - seekDRIVE))) / seekTAU #Seek Activity
        dsetp_dt = (-setp + F(Esetp * (0.9*nac + dlsSCALE * dls - setpDRIVE - setpDRIVE))) / setpTAU #Alcohol Variable
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - bingeDRIVE))) / bingeTAU #Binge Activity
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE))) / nacTAU #NAc Activity
        ddls_dt = (-dls + F(Edls * ( 3*dls + csTOdls * CS - dlsDRIVE)))/ dlsTAU
        dALCOHOL_dt = nac + dlsSCALE * dls # Alcohol consumed 
        dvta_dt = (-vta + F(Evta*( csTOvta * CS - vtaDRIVE))) / vtaTAU #VTA activity
       

        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, ddls_dt, dALCOHOL_dt, dvta_dt, dEnac_dt, dnacDRIVE_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    CS = np.heaviside(csDUR-t,0.5)
    y[4] = dlsSCALE * y[4]
    # y[8] = 0.9 * y[8]
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'CS':CS}

# Some extra functions defined for plotting 
def der_model(t, y , param): #This model is also defined separately, just for convenience. The system of ODEs is also defined in the "binge_model"
        seek, setp, binge, nac, dls,  ALCOHOL, vta, Enac, nacDRIVE = y

        csTOvta = param

        dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN - Enac)
        dnacDRIVE_dt = vta/nacdriveTAU + nacdriveDECAY*(driveMEAN - nacDRIVE)

        CS = np.heaviside(csDUR-t, 0.5) #Conditioned Stimulus
        dseek_dt = (-seek + F(Eseek * (binTOseek * binge + csTOseek * CS - spTOseek * setp - seekDRIVE))) / seekTAU #Seek Activity
        dsetp_dt = (-setp + F(Esetp * (0.9*nac + dlsSCALE * dls - setpDRIVE - setpDRIVE))) / setpTAU #Alcohol Variable
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - bingeDRIVE))) / bingeTAU #Binge Activity
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE))) / nacTAU #NAc Activity
        ddls_dt = (-dls + F(Edls * ( 3*dls + csTOdls * CS - dlsDRIVE)))/ dlsTAU
        dALCOHOL_dt = nac + dlsSCALE * dls # Alcohol consumed 
        dvta_dt = (-vta + F(Evta*( csTOvta * CS - vtaDRIVE))) / vtaTAU #VTA activity
       

        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, ddls_dt, dALCOHOL_dt, dvta_dt, dEnac_dt, dnacDRIVE_dt]


param_array = np.linspace(3, 5.5, 100)


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
    Enac = y['Int'][7]
    drive = y['Int'][8]
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
    axs[0,2].plot(t, nac + dls,  '--', label = 'Striatum', color = 'pink')
    axs[0,2].set_title('Striatal Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[0,2].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,2].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[0,2].set_ylim(0,1.1)
    axs[0,2].legend()

    # axs[1,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,0].plot(t, vta, label = 'DA', color = 'firebrick')
    axs[1,0].set_ylim(0,1)
    axs[1,0].set_title('VTA Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,0].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[1,0].legend()

    # axs[1,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,1].plot(t, CS, color = 'peru', label = 'CS')
    axs[1,1].plot(t, Enac, label = 'Enac' , color = 'blue')
    axs[1,1].plot(t, abs(drive) , color = 'pink', label = 'Nac Drive')

    # axs[1,1].plot(t, nacdrSCALE/Enac , color = 'orange', label = 'Nac Drive')

    # axs[1,1].plot(t, av, color = 'sienna', label = 'AV')
    axs[1,1].set_title('Conditioned Stimulus \n and NAc Parameters',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,1].set_xlabel('Time (min)',**tfont, fontsize='12')
    axs[1,1].legend()

    # axs[1,2].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,2].plot(t, alc, color = 'red')
    axs[1,2].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,2].set_ylabel('Volume (mL)',**tfont, fontsize='12')
    axs[1,2].set_xlabel('Time (min)',**tfont, fontsize='12')
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
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
        Excite = y['Int'][7]
        drive = y['Int'][8]
        CS = y['CS'] 
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
        axs[0,2].plot(t, nac +  dls,  '--', label = 'Striatum', color = 'pink')
        axs[0,2].set_title('Striatal Activity',**tfont, fontweight = 'bold', fontsize='14')
        axs[0,2].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[0,2].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[0,2].set_ylim(0,1.1)
        axs[0,2].legend()

        # axs[1,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
        axs[1,0].plot(t, vta, label = 'DA', color = 'firebrick')
        axs[1,0].set_ylim(0,1)
        axs[1,0].set_title('VTA Activity',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[1,0].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,0].legend()

        # axs[1,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
        axs[1,1].plot(t, CS, color = 'peru', label = 'CS')
        axs[1,1].plot(t, Excite, label = 'Enac' , color = 'blue')
        axs[1,1].plot(t, abs(drive) , color = 'orange', label = 'Nac Drive')
        # axs[1,1].plot(t, nacdrSCALE/Enac , color = 'pink', label = 'old Drive')


        # axs[1,1].plot(t, av, color = 'sienna', label = 'AV')
        axs[1,1].set_title('Conditioned Stimulus \n and NAc Parameters',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
        axs[1,1].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,1].legend()

        # axs[1,2].axvline(x=thresh, color = 'silver',linestyle='dashed')
        axs[1,2].plot(t, alc, color = 'red')
        axs[1,2].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
        axs[1,2].set_ylabel('Volume (mL)',**tfont, fontsize='12')
        axs[1,2].set_xlabel('Time (min)',**tfont, fontsize='12')
        axs[1,2].set_ylim(0,20)

        plt.subplots_adjust(wspace=0.3, hspace=0.4)
        return axs
   ani = animation.FuncAnimation(f, update, frames=len(param_array), interval=1,repeat=False)
   if save =='yes':
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulateNAcparam.gif', writer=writer)
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




# time = np.linspace(0,50,3)
time = np.array([0,15,30])
def td_vect(t,y0, time):
    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Seek', **tfont, fontsize = 15)
    ax.set_ylabel('Binge', **tfont, fontsize = 15)
    ax.set_zlabel('Time',rotation = 90, **tfont, fontsize = 15)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_zticks([0, 15/200, 30/200], ['0', '15', '30'], rotation=20, **tfont, fontsize = 12)

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    init = binge_model(time, y0, csTOvta)['Int'] # Solving the system with IC
    y_plot = binge_model(t, y0, csTOvta)['Int']
    new_traj = binge_model(t,[0.7, 0, 0.7, 0.1, 0, 0, 0, 0.6, -7],csTOvta)['Int']
    two_traj = binge_model(t,[0, 0, 1, 0.1, 0, 0, 0, 0.6, -7],csTOvta)['Int']

    seek_traj = new_traj[0]
    binge_traj = new_traj[2]

    seek_new = two_traj[0]
    binge_new = two_traj[2]

    seek_plot = y_plot[0]
    binge_plot = y_plot[2]
    for k in np.arange(len(time)):
        seek_array, binge_array , z = np.meshgrid(np.arange(0, 1.2, 0.2), np.arange(0, 1.2, 0.2), time[k]/200)
        dseek_array = np.zeros(seek_array.shape)
        dbinge_array = np.zeros(binge_array.shape)
        for i in np.arange(len(seek_array)):
            for j in np.arange(len(binge_array)):
                y = [init[0, k], init[1, k], init[2, k], init[3, k], init[4,k], init[5, k], init[6, k], init[7, k], init[8,k]]
                y[0] = seek_array[i, j,0]
                y[2] = binge_array[i, j,0]
                deriv = der_model(time, y, csTOvta)
                seek_d = deriv[0]
                dseek_array[i,j] = seek_d[k]
                dbinge_array[i,j] = deriv[2]
        ax.quiver(seek_array, binge_array, time[k]/200, dseek_array, dbinge_array, 0 ,length = 0.09, alpha = 0.6, arrow_length_ratio = 0.15)
        ax.plot3D(seek_plot, binge_plot, t/200, color = 'black', linewidth = 1.8)
        ax.plot3D(seek_traj, binge_traj, t/200, color = 'black', linewidth = 1.8)
        ax.plot3D(seek_new, binge_new, t/200, color = 'black', linewidth = 1.8)

        ax.scatter3D(seek_plot[0], binge_plot[0], t[0]/200, color = 'red', linewidth = 5)
        ax.scatter3D(seek_traj[0], binge_traj[0], t[0]/200, color = 'red', linewidth = 5)
        ax.scatter3D(seek_new[0], binge_new[0], t[0]/200, color = 'red', linewidth = 5)

        ax.scatter3D(seek_plot[-1], binge_plot[-1], t[-1]/200, color = 'black', linewidth = 2.5, marker = '^')
        # ax.scatter3D(seek_plot[int(len(t)/2)], binge_plot[int(len(t)/2)], t[int(len(t)/2)]/200, color = 'red', linewidth = 5)



    plt.show()
    return ax
sub_plots(t, y0, 'no', csTOvta)
# sub_plots_ani(t,y0, param_array, 'no')
td_vect(t,y0,time)