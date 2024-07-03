import numpy as np
from params import *
import matplotlib.pyplot as plt
import pythonizedModel 
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
tfont = {'fontname':'Times New Roman'}


#Plotting the Data
def sub_plots(t,y0):
    y= pythonizedModel.xppaut_model(t,y0)
    seek = y['Int'][0]
    binge = y['Int'][1]
    stop = y['Int'][2]
    nac = y['Int'][3]
    dls = y['Int'][4]
    alc = y['Int'][5]
    setp = 1-pythonizedModel.F(Esetp*(alc-TOLERANCE))
    vta = pythonizedModel.F(Evta*(alc - TOLERANCE*daFACTOR))
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
    axs[0,0].set_xlabel('Time (s)',**tfont, fontsize='12')
    axs[0,0].legend()

    axs[0,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[0,1].plot(t, binge, label ='Binge', color = 'mediumseagreen')
    axs[0,1].plot(t, stop, label = 'Stop', color = 'darkgreen')
    axs[0,1].set_title('Insular Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[0,1].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[0,1].set_xlabel('t (s)',**tfont, fontsize='12')
    axs[0,1].legend()

    axs[1,0].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,0].plot(t, vta, label = 'DA', color = 'lightcoral')
    axs[1,0].plot(t,nac, label = 'NAc', color = 'maroon')
    axs[1,0].plot(t,dls, label = 'DLS', color = 'crimson')
    axs[1,0].set_title('Subcortical Nuclei Activity',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,0].set_ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    axs[1,0].set_xlabel('t (s)',**tfont, fontsize='12')
    axs[1,0].legend()

    axs[1,1].axvline(x=thresh, color = 'silver',linestyle='dashed')
    axs[1,1].plot(t, alc, color = 'red')
    axs[1,1].set_title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='14')
    axs[1,1].set_ylabel('Volume (mL)',**tfont, fontsize='12')
    axs[1,1].set_xlabel('t (s)',**tfont, fontsize='12')
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    plt.show()

def comb_plots(t,y0):
    y= pythonizedModel.xppaut_model(t,y0)
    seek = y['Int'][0]
    binge = y['Int'][1]
    stop = y['Int'][2]
    nac = y['Int'][3]
    dls = y['Int'][4]
    alc = y['Int'][5]
    setp = 1-pythonizedModel.F(Esetp*(alc-TOLERANCE))
    vta = pythonizedModel.F(Evta*(alc - TOLERANCE*daFACTOR))
    for n in np.arange(len(alc)):
        if alc[n]>=TOLERANCE:
            thresh = t[n] #time at which threshold is reached 
            index = n #index of when threshold is reached
            break
    plt.figure(figsize=(15,5))
    plt.plot(t, seek, label = 'Seek', color = 'midnightblue', linewidth='1.2')
    plt.plot(t,setp, label = 'Setpoint', color = 'royalblue', linewidth='1.2')
    plt.plot(t, binge, label ='Binge', color = 'mediumseagreen', linewidth='1.2')
    plt.plot(t, stop, label = 'Stop', color = 'darkgreen', linewidth='1.2')
    plt.plot(t, vta, label = 'DA', color = 'lightcoral', linewidth='1.2')
    plt.plot(t,nac, label = 'NAc', color = 'maroon', linewidth='1.2')
    plt.plot(t,dls, label = 'DLS', color = 'crimson', linewidth='1.2')
    plt.xlabel('t (s)',**tfont, fontsize='15',fontweight = 'bold')
    plt.ylabel('Firing Rate (Hz)',**tfont, fontsize='15',fontweight = 'bold')
    plt.legend()
    plt.title('Simplified Model of Front-Loading',**tfont, fontweight = 'bold', fontsize='20')
    plt.show()

def ind_plots(t,y0):
    y= pythonizedModel.xppaut_model(t,y0)['Int']
    seek = y[0]
    binge = y[1]
    stop = y[2]
    nac = y[3]
    dls = y[4]
    alc = y[5]
    setp = 1-pythonizedModel.F(Esetp*(alc-TOLERANCE))
    vta = pythonizedModel.F(Evta*(alc - TOLERANCE*daFACTOR))
    for n in np.arange(len(alc)):
        if alc[n]>=TOLERANCE:
            thresh = t[n] #time at which threshold is reached 
            index = n #index of when threshold is reached
            break
    
    plt.axvline(x=thresh, color = 'silver',linestyle='dashed')
    plt.plot(t, seek+setp, label = 'Combined', color = 'lightsteelblue')
    plt.plot(t, seek, label = 'Seek', color = 'midnightblue')
    plt.plot(t, setp, label = 'Setpoint', color = 'royalblue')
    plt.title('mPFC Activity', **tfont, fontweight = 'bold', fontsize='14')
    plt.ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    plt.xlabel('t (s)',**tfont, fontsize='12')
    plt.legend()
    plt.show()

    plt.axvline(x=thresh, color = 'silver',linestyle='dashed')
    plt.plot(t, binge, label ='Binge', color = 'mediumseagreen')
    plt.plot(t, stop, label = 'Stop', color = 'darkgreen')
    plt.title('Insular Activity',**tfont, fontweight = 'bold', fontsize='14')
    plt.ylabel('Firing Rate (Hz)',**tfont, fontsize='12')
    plt.xlabel('t (s)',**tfont, fontsize='12')
    plt.legend()
    plt.show()
      
    plt.axvline(x=thresh, color = 'silver',linestyle='dashed')
    plt.plot(t, vta, label = 'DA', color = 'lightcoral')
    plt.plot(t,nac, label = 'NAc', color = 'maroon')
    plt.plot(t,dls, label = 'DLS', color = 'crimson')
    plt.title('Subcortical Nuclei Activity',**tfont, fontweight = 'bold', fontsize='18')
    plt.ylabel('Firing Rate (Hz)',**tfont, fontsize='15')
    plt.xlabel('t (s)',**tfont, fontsize='15')
    plt.legend()
    plt.show()
    
    plt.axvline(x=thresh, color = 'silver',linestyle='dashed')
    plt.plot(t, alc, color = 'red')
    plt.title('Alcohol Consumed',**tfont, fontweight = 'bold', fontsize='18')
    plt.ylabel('Volume (mL)',**tfont, fontsize='15')
    plt.xlabel('t (s)',**tfont, fontsize='15')
    plt.show()
    plt.show()


#Plotting the Nullcline and Phase-Space

def phase_space(y0, y_traj, t,n,m,name): 
    #n, m = index for neuron group for phase plane analysis 
    #y0 = [seek, binge, stop, nac, dls, alc]
    #name = ['group one', 'group two']
    
    #Defining the shape of the grid
    x1min=-0.1
    x1max=1
    numptsx1=12 
    x1list=np.linspace(x1min,x1max,numptsx1)
    x2min=-0.1
    x2max=1
    numptsx2=12 
    x2list=np.linspace(x2min,x2max,numptsx2)
    x1array,x2array = np.meshgrid(x1list,x2list)
    dx1dt_array=np.zeros(x1array.shape)
    dx2dt_array=np.zeros(x1array.shape)
    
    #Calculating and Plotting the Vector Fields 
    for i in np.arange(numptsx1):
        for j in np.arange(numptsx2):      
            y0[n] = x1array[i,j]
            y0[m] = x2array[i,j]
            deriv = pythonizedModel.derModel(0,y0)
            dx1dt_array[i,j]= deriv[n]
            dx2dt_array[i,j]= deriv[m]
    fig = plt.figure(figsize=(8, 6))
    plt.quiver(x1array,x2array,dx1dt_array,dx2dt_array,alpha=0.3) # quiver: plots vector field

    #Calculating the Nullclines
    x1list_fine=np.linspace(x1min,x1max,250)
    x2list_fine=np.linspace(x2min,x2max,250)
    nullcline = np.zeros((5,250))
    for i in np.arange(250):
        y0[n] = x1list_fine[i]
        y0[m] = x2list_fine[i]
        setp = 1-pythonizedModel.F(Esetp * (y0[4] - TOLERANCE))
        vta = pythonizedModel.F(Evta*(y0[5] - TOLERANCE*daFACTOR))
        aps = Eaps*(y0[3] + y0[4])/2
        nullcline[0,i] = (pythonizedModel.F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
        nullcline[1,i] = (pythonizedModel.F(Ebinge * (stopTObin * y0[2] - seekTObin * y0[0] - bingeDRIVE))) / bingeTAU
        nullcline[2,i] = (pythonizedModel.F(Estop * (binTOstop * y0[1] - spTOstop * setp - stopDRIVE))) / stopTAU
        nullcline[3,i] = (pythonizedModel.F(Enac * (-vtaTOnac * vta - seekTOnac * y0[0] - binTOnac * y0[1] - nacDRIVE))) / nacTAU
        nullcline[4,i] = (pythonizedModel.F(Edls * (-binTOdls * y0[1] - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
    
    plt.plot(x1list_fine,nullcline[n,:],'b-',alpha=0.8,linewidth = 2.5)
    plt.plot(nullcline[m,:],x2list_fine,'b-',alpha=0.8,linewidth = 2.5)
    plt.xlabel(name[0], **tfont, fontweight = 'bold', fontsize = '15')
    plt.ylabel(name[1], **tfont, fontweight = 'bold', fontsize = '15')

    #Plotting the Trajectories
    t_span = (0, 50)
    y= pythonizedModel.xppaut_model(t,y_traj)['Int']
    traj = np.zeros((6,len(t)))
    traj[0,:] = y[0] #seek
    traj[1,:] = y[1] #binge
    traj[2,:] = y[2] #stop
    traj[3,:] = y[3] #nac
    traj[4,:] = y[4] #dls
    traj[5,:] = y[5] #alc
    setp = 1-pythonizedModel.F(Esetp*(y[5]-TOLERANCE))
    vta = pythonizedModel.F(Evta*(y[5] - TOLERANCE*daFACTOR))
    for i in np.arange(len(y[5])):
        if y[5][i]>=TOLERANCE:
            thresh = t[i] #time at which threshold is reached 
            index = i #index of when threshold is reached
            break
    print(len(traj[1,:]))
    print(len(traj[2,:]))
    print('n=',n)
    plt.plot(traj[n,:], traj[m,:], color = 'black')
    plt.plot(traj[n,0], traj[m,0], markersize = '15' , marker = '*', color = 'lightcoral')
    plt.plot(traj[n,index-10], traj[m,index-10], markersize = '10' , marker = 'o', color = 'lightcoral')
    plt.plot(traj[n,-1], traj[m,-1], markersize = '10' , marker = 'o', color = 'red')
    plt.show()
