import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from mpl_toolkits import mplot3d
tfont = {'fontname':'Times New Roman'}

Ebinge = 2
Esetp = 1.84
Eseek= 2
Evta = 6
Edls = 2

seekTAU = 1
bingeTAU = .5
nacTAU = 1
setpTAU = 20
vtaTAU = 1
dlsTAU = 4

seekDRIVE = 0
bingeDRIVE = -1
setpDRIVE = -1
vtaDRIVE = -1.4
dlsDRIVE = -2

spTOseek = 13
seekTOnac = 4
seekTObin = 2.5
binTOseek = 2.5
binTOnac = 1
vtaTOnac = 2.5
csTOseek = 4
csTOvta = 2.5 # modulate this connection to change magnitude of DA peak (1 to 3)
csTOdls = 3
nacTOsetp = .1 #length of front-loading
dlsTOdls = 5


csDUR = 3
EnacDECAY = 0.001
nacDriveDECAY = 0.003
EnacTAU = 1
nacDriveTAU = 0.7
EnacMEAN = 0.6
driveMEAN = -7
dlsSCALE = 0.1

y0 = [0.1, 0, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]
t = np.linspace(0,120,500)
#setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

def xppaut_model(t,y, csTOvta=csTOvta):
    def model(t, y):
        setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE = y
        cs = np.heaviside(csDUR-t, .5)

        dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN-Enac)
        dnacDRIVE_dt =  vta/nacDriveTAU + nacDriveDECAY*(driveMEAN-nacDRIVE)
        
        dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsSCALE*dls) + setpDRIVE))) / setpTAU
        dseek_dt = (-seek + F(Eseek * (-spTOseek * setp + csTOseek * cs + binTOseek * binge + seekDRIVE))) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek + bingeDRIVE))) / bingeTAU
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE))) / nacTAU
        ddls_dt = (-dls + F(Edls * (dlsTOdls * dls + csTOdls * cs + dlsDRIVE))) / dlsTAU
        
        dvta_dt = (-vta + F(Evta*(csTOvta * cs + vtaDRIVE))) / vtaTAU
        dALCOHOL_dt = (nac+dlsSCALE*dls)
        return [dsetp_dt, dseek_dt, dbinge_dt, dnac_dt, ddls_dt, dvta_dt, dALCOHOL_dt, dEnac_dt, dnacDRIVE_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    y[4] = dlsSCALE*y[4]
    cs = np.heaviside(csDUR-t, .5)
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'CS':cs}
def der_model(t, y):
    setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE = y
    cs = np.heaviside(csDUR-t, .5)

    dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN-Enac)
    dnacDRIVE_dt =  vta/nacDriveTAU + nacDriveDECAY*(driveMEAN-nacDRIVE)
    
    dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsSCALE*dls) + setpDRIVE))) / setpTAU
    dseek_dt = (-seek + F(Eseek * (-spTOseek * setp + csTOseek * cs + binTOseek * binge + seekDRIVE))) / seekTAU
    dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek + bingeDRIVE))) / bingeTAU
    dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE))) / nacTAU
    ddls_dt = (-dls + F(Edls * (dlsTOdls * dls + csTOdls * cs + dlsDRIVE))) / dlsTAU
    
    dvta_dt = (-vta + F(Evta*(csTOvta * cs + vtaDRIVE))) / vtaTAU
    dALCOHOL_dt = (nac+dlsSCALE*dls)
    return [dsetp_dt, dseek_dt, dbinge_dt, dnac_dt, ddls_dt, dvta_dt, dALCOHOL_dt, dEnac_dt, dnacDRIVE_dt]

def runGraphs(time=120, fuzz=False, save=False, anim=False):
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    t = np.linspace(0, time, 300)
    y = xppaut_model(t, fuzz)
    sp = axs[0, 0].plot(t, y['Int'][0], label="Setpoint", color = 'royalblue')[0]
    seek = axs[0, 0].plot(t, y['Int'][1], label="Seek", color = 'midnightblue')[0]
    comb = axs[0, 0].plot(t, (y['Int'][0]+y['Int'][1])/2, '--', label="mPFC Average", color = 'lightblue')[0]
    axs[0,0].set_title("mPFC Activity")
    
    bin = axs[0, 1].plot(t, y['Int'][2], label="Binge", color = 'mediumseagreen')[0]
    axs[0,1].set_title("Insular Activity")
    
    nac = axs[0, 2].plot(t, y['Int'][3], label="NAc", color = 'maroon')[0]
    dls = axs[0, 2].plot(t, y['Int'][4], label="DLS", color='red')[0]
    sum = axs[0,2].plot(t, y['Int'][3]+y['Int'][4], '--',label="Striatum", color='tomato')[0]
    axs[0, 2].set_title("Striatal Activity")

    da = axs[1, 0].plot(t, y['Int'][5], label="VTA", color = 'lightcoral')[0]
    axs[1,0].set_title("VTA Activity")

    excNac = axs[1,1].plot(t, y['Int'][7], label="$E_{NAc}$")[0]
    driNac = axs[1, 1].plot(t, abs(y['Int'][8]), label='$drive_{NAc}$')[0]
    axs[1,1].set_title("DA Modulation of NAc Parameters")

    
    alc = axs[1, 2].plot(t, y['Int'][6], label='Alcohol Vol.', color = 'red')[0]
    axs[1,2].set_title("Alcohol Consumption")
    for z in range(3):
         axs[0,z].set_ylabel('Firing Rate (Hz)')
    axs[1,0].set_ylabel('Firing Rate (Hz)')
    axs[1,2].set_ylabel('Volume')
    for i in range(2):
         for j in range(3):
            axs[i, j].legend()
            axs[i,j].set_xlabel('Time (mins)')
            axs[i,j].fill_between(t, 30, where=[(t >= 0) and (t <= 3) for t in t], color = 'grey', alpha = 0.15, linewidth = 0.05) 
            axs[i,j].set_xlim(0,t[-1])
    axs[0, 0].set_ylim(0, 1)
    axs[0, 1].set_ylim(0, 1)
    axs[0, 2].set_ylim(0, 1.2)
    axs[1, 0].set_ylim(0, 1)
    axs[1, 1].set_ylim(0, 7.5)
    axs[1, 2].set_ylim(0, 30)

    frames=100
    def update(frame):
        y = xppaut_model(t, fuzz, csTOvta=2*(frame/frames)+0)
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
        driNac.set_data(t, y['Int'][8])

        return (sp, seek, nac)
    if anim:
        ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=10, repeat = False)
    plt.tight_layout()
    if save:
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulation.gif', writer=writer)
    plt.show()

def td_vect(t,y0, time):
    time = np.array([0,20,40])
    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Seek', **tfont, fontsize = 15)
    ax.set_ylabel('Binge', **tfont, fontsize = 15)
    ax.set_zlabel('Time',rotation = 90, **tfont, fontsize = 15)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_zticks([0, 15/200, 30/200], ['0', '20', '40'], rotation=20, **tfont, fontsize = 12)

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    t = np.linspace(0,time[-1], 500)
    init = xppaut_model(time, y0, csTOvta)['Int'] # Solving the system with IC
    y_plot = xppaut_model(t, y0, csTOvta)['Int']
    new_traj = xppaut_model(t,[0.1, 0.9, 0.9, 0.1, 0, 0, 0, EnacMEAN, driveMEAN],csTOvta)['Int']
    two_traj = xppaut_model(t,[0.1, 0.2, 0.8, 0.1, 0, 0, 0, EnacMEAN, driveMEAN],csTOvta)['Int']

    seek_traj = xppaut_model(t,[0.1, 0.6, 0.6, 0.1, 0, 0, 0, EnacMEAN, driveMEAN],csTOvta)['Int'][1]
    binge_traj = xppaut_model(t,[0.1, 0.6, 0.6, 0.1, 0, 0, 0, EnacMEAN, driveMEAN],csTOvta)['Int'][2]

    seek_new = two_traj[1]
    binge_new = two_traj[2]

    seek_plot = y_plot[1]
    binge_plot = y_plot[2]
    for k in np.arange(len(time)):
        seek_array, binge_array , z = np.meshgrid(np.arange(0, 1.2, 0.2), np.arange(0, 1.2, 0.2), time[k]/200)
        dseek_array = np.zeros(seek_array.shape)
        dbinge_array = np.zeros(binge_array.shape)
        for i in np.arange(len(seek_array)):
            for j in np.arange(len(binge_array)):
                y = [init[0, k], init[1, k], init[2, k], init[3, k], init[4,k], init[5, k], init[6, k], init[7, k], init[8,k]]
                y[1] = seek_array[i, j, 0]
                y[2] = binge_array[i, j, 0]
                deriv = der_model(time, y)
                seek_d = deriv[1]
                dseek_array[i,j] = seek_d[k]
                dbinge_array[i,j] = deriv[2]
        # ax.scatter3D(seek_plot[0], binge_plot[0], t[0]/200, color = 'red', linewidth = 5)
        ax.scatter3D(seek_traj[0], binge_traj[0], t[0]/200, color = 'red', linewidth = 5)
        # ax.scatter3D(seek_new[0], binge_new[0], t[0]/200, color = 'red', linewidth = 5)
        ax.quiver(seek_array, binge_array, time[k]/200, dseek_array, dbinge_array, 0 ,length = 0.1, alpha = 0.6, arrow_length_ratio = 0.15)
        # ax.plot3D(seek_plot, binge_plot, t/200, color = 'black', linewidth = 1.8)
        # ax.plot3D(seek_traj, binge_traj, t/2500, color = 'black', linewidth = 1.8)
        # ax.plot3D(seek_new, binge_new, t/200, color = 'black', linewidth = 1.8)    

        # ax.scatter3D(seek_plot[-1], binge_plot[-1], t[-1]/200, color = 'black', linewidth = 2.5, marker = '^')
    plt.show()
    return ax

def ind_plots(graph, t=t):
    #graph: PFC, Insula, STR, VTA, Alc, Param
    y = xppaut_model(t, y0)['Int']
    f, ax = plt.subplots(1, 2, sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})
    
    if graph == 'PFC':
        for z in range(2):
            ax[z].plot(t,y[0], label="Setpoint", color = 'royalblue')
            ax[z].plot(t, y[1], label="Seek", color = 'midnightblue')
            ax[z].plot(t, (y[0]+y[1])/2, '--', label="mPFC Average", color = 'lightblue')
            f.suptitle('mPFC Activity', fontsize = 15, fontweight = 'bold')
    if graph == 'Insula':
        for z in range(2):
            ax[z].plot(t, y[2], label="Binge", color = 'mediumseagreen')
            f.suptitle('Insular Activity', fontsize = 15, fontweight = 'bold')
    if graph == 'STR':
        for z in range(2):
            ax[z].plot(t, y[3], label="NAc", color = 'maroon')
            ax[z].plot(t, y[4], label="DLS", color='red')
            ax[z].plot(t, y[3]+y[4], '--',label="Striatum", color='tomato')
            f.suptitle('Striatal Activity', fontsize = 15, fontweight = 'bold')
    if graph == 'VTA':
        for z in range(2):
            ax[z].plot(t, y[5], label="VTA", color = 'lightcoral')
            f.suptitle('VTA Activity', fontsize = 15, fontweight = 'bold')
    if graph == 'Alc':
        f = plt.figure()
        plt.plot(t, y[6], label='Alcohol Vol.', color = 'red')
        plt.xlabel('Time (mins)')
        plt.ylabel('Volume')
        plt.title('Alcohol Consumption', fontsize = 15, fontweight = 'bold')
        plt.legend()
        plt.show()
    if graph == "Param":
        f = plt.figure()
        plt.plot(t, y[7], label="$E_{NAc}$")[0]
        plt.plot(t, abs(y[8]), label='|$drive_{NAc}$|')[0]
        plt.title("DA Modulation of NAc Parameters", fontsize = 15, fontweight = 'bold')
        plt.legend()
        plt.xlabel('Time (mins)')
        plt.show()

    if graph == 'PFC' or 'Insula' or 'STR' or 'VTA':
        ax[0].set_xlim(0, 40)
        ax[1].set_xlim(115, 120)
        ax[0].spines['right'].set_visible(False)
        ax[1].spines['left'].set_visible(False)
        ax[1].tick_params(left = False) 
        d = .015  
        kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
        ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
        ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
        ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
        ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
        ax[0].set_xlabel('Time (mins)')
        ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10)
        ax[0].xaxis.set_label_coords(0.6, -0.09)
        ax[1].legend()
        plt.show()

def DA_graphs(F, S, M, L):
    fail = xppaut_model(t,y0,F)['Int']
    low = xppaut_model(t,y0,S)['Int']
    medium = xppaut_model(t,y0,M)['Int']
    high = xppaut_model(t,y0, L)['Int']
    
    fig, ax = plt.subplots(1,2, figsize = (10,5))
    ax[0].plot(t, fail[5], label = 'No DA')
    ax[0].plot(t, low[5], label = 'Low DA')
    ax[0].plot(t, medium[5], label = 'Medium DA')
    ax[0].plot(t, high[5], label = 'High DA')
    ax[0].set_xlim(0,10)
    ax[0].set_xlabel('Time (mins)', fontsize ='12')
    ax[0].set_ylabel('Firing Rate (Hz)', fontsize ='12') 
    ax[1].plot(t, fail[6], label = 'No DA')
    ax[1].plot(t, low[6], label = 'Low DA')
    ax[1].plot(t, medium[6], label = 'Medium DA')
    ax[1].plot(t, high[6], label = 'High DA')
    ax[1].set_ylim(0,30)
    ax[1].set_xlabel('Time (mins)', fontsize ='12')
    ax[1].set_ylabel('Volume', fontsize ='12')
    ax[1].legend()
    plt.tight_layout()
    plt.show()
    
    
    param_array = np.linspace(0, 2.5, 100)
    final_alcohol = []
    peak_vta = []
    peak_nac = []

    for n in np.arange(len(param_array)):
         y= xppaut_model(t, y0, param_array[n])['Int']
         alc = y[6]
         final_alcohol.append(alc[-1])
         peak_vta.append(np.max(y[5]))
         peak_nac.append(np.max(y[3]))

    fig, ax1 = plt.subplots(figsize = (8,6))
    color = 'tab:red'
    ax1.set_xlabel('Peak VTA Activity (Hz)', fontsize = '12')
    ax1.set_ylabel('Total Alcohol Consumed', color=color, fontsize = '12')
    ax1.scatter(peak_vta, final_alcohol, color=color, linewidths = 0.05)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  
    color = 'tab:blue'
    ax2.set_ylabel('Peak NAc Activity (Hz)', color=color, fontsize = '12', rotation = -90, labelpad = 20)  
    ax2.scatter(peak_vta, peak_nac, color=color, linewidths = 0.05)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('Effect of DA Release', fontsize = '15', fontweight = 'bold')
    fig.tight_layout()  
    plt.show()


DA_graphs(0.5, 1.5, 1.65, 2)
ind_plots('PFC') #graph: PFC, Insula, STR, VTA, Alc, Param
# td_vect(t, y0)
# runGraphs(120, anim=True)









# low = xppaut_model(t,y0, 1)['Int'][6]
# low_vta = xppaut_model(t,y0, 1)['Int'][5]
# medium = xppaut_model(t,y0,1.65)['Int'][6]
# medium_vta = xppaut_model(t,y0,1.65)['Int'][5]
# high = xppaut_model(t,y0, 2)['Int'][6]
# high_vta = xppaut_model(t,y0, 2)['Int'][5]

# fig, ax = plt.subplots(1,2, figsize = (8,5))
# ax[1].plot(t, low, label = 'Low DA')
# ax[1].plot(t, medium, label = 'Medium DA')
# ax[1].plot(t, high, label = 'High DA')
# ax[1].set_xlabel('Time (mins)')
# ax[1].set_ylabel('Volume')
# ax[1].legend()
# ax[0].plot(t, low_vta, label = 'Low DA')
# ax[0].plot(t, medium_vta, label = 'Medium DA')
# ax[0].plot(t, high_vta, label = 'High DA')
# ax[0].set_xlim(0,10)
# ax[0].set_xlabel('Time (mins)')
# ax[0].set_ylabel('VTA')

# plt.show()






# param_array = np.linspace(0, 2.5, 100)
# final_alcohol = []
# peak_vta = []
# peak_nac = []

# for n in np.arange(len(param_array)):
#      y= xppaut_model(t, y0, param_array[n])['Int']
#      alc = y[6]
#      final_alcohol.append(alc[-1])
#      peak_vta.append(np.max(y[5]))
#      peak_nac.append(np.max(y[3]))

# fig, ax1 = plt.subplots(figsize = (8,6))
# color = 'tab:red'
# ax1.set_xlabel('Peak VTA Activity (Hz)', fontsize = '12')
# ax1.set_ylabel('Total Alcohol Consumed', color=color, fontsize = '12')
# ax1.scatter(peak_vta, final_alcohol, color=color, linewidths = 0.05)
# ax1.tick_params(axis='y', labelcolor=color)

# ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis
# color = 'tab:blue'
# ax2.set_ylabel('Peak NAc Activity (Hz)', color=color, fontsize = '12', rotation = -90, labelpad = 20)  # we already handled the x-label with ax1
# ax2.scatter(peak_vta, peak_nac, color=color, linewidths = 0.05)
# ax2.tick_params(axis='y', labelcolor=color)

# plt.title('Effect of DA Release', fontsize = '15', fontweight = 'bold')
# fig.tight_layout()  # otherwise the right y-label is slightly clipped
# plt.show()


