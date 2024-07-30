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
ytraj = [0.1, 0, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]

t = np.linspace(0,120,500)
#setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

x= np.linspace(-10,10,500)
def sig_graphs(x):
    fig = plt.figure(figsize = (9,8))
    plt.plot(x, F(x), 'k', linewidth = 7)
    plt.plot(x,F(-x), '--', color = 'black', linewidth = 7)
    plt.xlabel('Inputs', fontsize = 28)
    plt.ylabel('F', fontsize = 28)
    plt.xticks(np.array([-10,-5,0,5,10]),fontsize = 25)
    plt.yticks(np.array([0,0.5,1]),fontsize = 25)
    plt.show()
    fig.savefig('/Users/amyrude/Downloads/sig.png', transparent=True, dpi=350)

# sig_graphs(x)


def xppaut_model(t,y0, noise = False, csTOvta=csTOvta, seekTObin=seekTObin, binTOseek=binTOseek):
    def model(t, y0):
        setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE = y0
        if noise ==True:
            cs = np.heaviside(csDUR-t, .5)

            dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN-Enac)
            dnacDRIVE_dt =  vta/nacDriveTAU + nacDriveDECAY*(driveMEAN-nacDRIVE)
            
            dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsSCALE*dls) + setpDRIVE)+ np.random.normal(0,0.3)) ) / setpTAU
            dseek_dt = (-seek + F(Eseek * (-spTOseek * setp + csTOseek * cs + binTOseek * binge + seekDRIVE)+ np.random.normal(0,0.3))) / seekTAU
            dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek + bingeDRIVE)+ np.random.normal(0,0.3))) / bingeTAU
            dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE)+ np.random.normal(0,0.3))) / nacTAU
            ddls_dt = (-dls + F(Edls * (dlsTOdls * dls + csTOdls * cs + dlsDRIVE)+ np.random.normal(0,0.3))) / dlsTAU
            
            dvta_dt = (-vta + F(Evta*(csTOvta * cs + vtaDRIVE)+ np.random.normal(0,0.3))) / vtaTAU
            dALCOHOL_dt = (nac+dlsSCALE*dls)
           
        else:
            cs = np.heaviside(csDUR-t, .5)

            dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN-Enac)
            dnacDRIVE_dt =  vta/nacDriveTAU + nacDriveDECAY*(driveMEAN-nacDRIVE)
            
            dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsSCALE*dls) + setpDRIVE)) ) / setpTAU
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


def der_model(t, y, seekTObin=seekTObin, binTOseek=binTOseek):
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



def runGraphs(time=120, fuzz=False, save_sub=False, save_ani = False, anim=False):
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    t = np.linspace(0, time, 300)
    y = xppaut_model(t, fuzz)
    sp = axs[0, 0].plot(t, y['Int'][0], label="Setpoint", color = 'royalblue')[0]
    seek = axs[0, 0].plot(t, y['Int'][1], label="Seek", color = 'midnightblue')[0]
    comb = axs[0, 0].plot(t, (y['Int'][0]+y['Int'][1])/2, '--', label="mPFC Average", color = 'lightblue')[0]
    axs[0,0].set_title("mPFC Activity", fontsize = 15, fontweight = 'bold')
    
    bin = axs[0, 1].plot(t, y['Int'][2], label="Binge", color = 'mediumseagreen')[0]
    axs[0,1].set_title("Insular Activity", fontsize = 15, fontweight = 'bold')
    
    nac = axs[0, 2].plot(t, y['Int'][3], label="NAc", color = 'maroon')[0]
    dls = axs[0, 2].plot(t, y['Int'][4], label="DLS", color='red')[0]
    sum = axs[0,2].plot(t, y['Int'][3]+y['Int'][4], '--',label="Striatum", color='tomato')[0]
    axs[0, 2].set_title("Striatal Activity", fontsize = 15, fontweight = 'bold')

    da = axs[1, 0].plot(t, y['Int'][5], label="VTA", color = 'lightcoral')[0]
    axs[1,0].set_title("VTA Activity", fontsize = 15, fontweight = 'bold')

    excNac = axs[1,1].plot(t, y['Int'][7], label="$E_{NAc}$")[0]
    driNac = axs[1, 1].plot(t, abs(y['Int'][8]), label='$|drive_{NAc}|$')[0]
    axs[1,1].set_title("DA Modulation of \n NAc Parameters", fontsize = 15, fontweight = 'bold')

    
    alc = axs[1, 2].plot(t, y['Int'][6], label='Alcohol Vol.', color = 'red')[0]
    axs[1,2].set_title("Alcohol Consumption", fontsize = 15, fontweight = 'bold')
    for z in range(3):
         axs[0,z].set_ylabel('Firing Rate (Hz)', fontsize = 13)
    axs[1,0].set_ylabel('Firing Rate (Hz)', fontsize = 13)
    axs[1,2].set_ylabel('Volume', fontsize = 13)
    for i in range(2):
         for j in range(3):
            axs[i,j].set_xlabel('Time (mins)', fontsize = 13)
            axs[i,j].fill_between(t, 30, where=[(t >= 0) and (t <= 3) for t in t], color = 'grey', alpha = 0.15, linewidth = 0.05, label = 'CS') 
            axs[i,j].set_xlim(0,120)
            axs[i, j].legend()

    axs[0, 0].set_ylim(0, 1)
    axs[0, 1].set_ylim(0, 1)
    axs[0, 2].set_ylim(0, 1.2)
    axs[1, 0].set_ylim(0, 1)
    axs[1, 1].set_ylim(0, 7.5)
    axs[1, 2].set_ylim(0, 30)
    
    frames=300
    def update(frame):
        y = xppaut_model(t, fuzz, csTOvta=3.5*(frame/frames)+0) # seekTObin=1*(frame/frames)+2.6, binTOseek=1*(frame/frames)+2.6
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
        driNac.set_data(t, abs(y['Int'][8]))

        return (sp, seek, nac)
    if anim:
        ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=20, repeat = False)
    plt.tight_layout()
    # plt.subplots_adjust(wspace=0.3, hspace=0.4)
    if save_sub == True:      
        fig.savefig('/Users/amyrude/Downloads/subplot.png', transparent=True, dpi=350)
    if save_ani==True:
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulation.gif', writer=writer)
    plt.show()



def td_vect(t,y0):
    time = np.array([0,20,40])
    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Seek',  fontsize = 20, labelpad = 10)
    ax.set_ylabel('Binge',  fontsize = 20, labelpad = 13)
    ax.set_zlabel('Time',rotation = 90, fontsize = 20)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize = 16)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ['0', '0.25', '0.5', '0.75', '1.0'],   fontsize = 16)
    ax.set_zticks([0, 20/200, 40/200], ['0', '20', '40'], rotation=20,  fontsize = 16)

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    t = np.linspace(0,time[-1], 500)
    init = xppaut_model(time, y0, csTOvta)['Int'] # Solving the system with IC

    for k in np.arange(len(time)):
        seek_array, binge_array , z = np.meshgrid(np.arange(0, 1.2, 0.2), np.arange(0, 1.2, 0.2), time[k]/200)
        dseek_array = np.zeros(seek_array.shape)
        dbinge_array = np.zeros(binge_array.shape)
        for i in np.arange(len(seek_array)):
            for j in np.arange(len(binge_array)):
                y = [init[0, k], init[1, k], init[2, k], init[3, k], init[4,k], init[5, k], init[6, k], init[7, k], init[8,k]]
                y[1] = seek_array[i, j, 0]
                y[2] = binge_array[i, j, 0]
                deriv = der_model(time[k],y)
                dseek_array[i,j] = deriv[1]
                dbinge_array[i,j] =  deriv[2]
        ax.quiver(seek_array, binge_array, time[k]/200, dseek_array, dbinge_array, 0 ,length = 0.09, alpha = 0.75, arrow_length_ratio = 0.25)
 
    
    traj_one = xppaut_model(t, y0)['Int']
    traj_two = xppaut_model(t, np.array([0.1, 0.2, 0.6, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]))['Int']
    traj_three = xppaut_model(t, [0.1, 0.6, 0.2, 0.1, 0, 0, 0, EnacMEAN, driveMEAN])['Int']

    
    ax.scatter3D(traj_one[1][0], traj_one[2][0], t[0]/200, color = 'red', linewidth = 5.5)
    ax.plot3D(traj_one[1], traj_one[2], t/200, color = 'black', linewidth =2)
    ax.scatter3D(traj_two[1][0], traj_two[2][0], t[0]/200, color = 'red', linewidth = 5.5)
    ax.plot3D(traj_two[1], traj_two[2], t/200, color = 'black', linewidth = 2)
    ax.scatter3D(traj_three[1][0], traj_three[2][0], t[0]/200, color = 'red', linewidth = 5.5)
    ax.plot3D(traj_three[1], traj_three[2], t/200, color = 'black', linewidth = 2)
    ax.scatter3D(traj_one[1][-1], traj_one[2][-1], t[-1]/200, color = 'black', linewidth = 2.5, marker = '^')
    ax.view_init(elev=13., azim=-109)
    fig.savefig('/Users/amyrude/Downloads/vect_field.png', transparent=True, dpi=350)
    plt.show()

    return ax


y = xppaut_model(t, y0, noise=True)['Int']
def ind_plots(graph, csTOvta=csTOvta, t=t, save=True):
    #graph: PFC, Insula, STR, VTA, Alc, Param
    f, ax = plt.subplots(1, 2, figsize = (7,7), sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})

    if graph == 'PFC':

        for z in range(2):
            ax[z].plot(t,y[0], label="Setpoint", color = 'royalblue', linewidth = 4)
            ax[z].plot(t, y[1], label="Seek", color = 'midnightblue', linewidth = 4)
            ax[z].plot(t, (y[0]+y[1])/2, '--', label="mPFC Average", color = 'lightblue', linewidth = 4)
            f.suptitle('mPFC Activity', fontsize = 20, fontweight = 'bold')
            
    if graph == 'Insula':
        for z in range(2):
            ax[z].plot(t, y[2], label="Binge", color = 'mediumseagreen')
            f.suptitle('Insular Activity', fontsize = 15, fontweight = 'bold')
    if graph == 'Cortex':
        for z in range(2):
            ax[z].plot(t, y[2], label="Binge", color = 'darkgreen', linewidth = 4)
            ax[z].plot(t,y[0], label="Setpoint", color = 'royalblue', linewidth = 4)
            ax[z].plot(t, y[1], label="Seek", color = 'midnightblue', linewidth = 4)
            ax[z].plot(t, (y[0]+y[1])/2, '--', label="mPFC \n Average", color = 'lightblue', linewidth = 4)
            # f.suptitle('Cortical Activity', fontsize = 20, fontweight = 'bold')

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
    if graph == 'Subcort':
        for z in range(2):
            ax[z].plot(t, y[3], label="NAc", color = 'chocolate', linewidth = 4)
            ax[z].plot(t, y[4], label="DLS", color='darkorange', linewidth = 4)
            ax[z].plot(t, y[3]+y[4], '--',label="Striatum", color='tan', linewidth = 4)
            ax[z].plot(t, y[5], label="VTA", color = 'indigo', linewidth = 4)


    if graph == 'Alc':
        f = plt.figure(figsize =(7,7)) #figsize = (14,6)
        plt.plot(t, y[6],  color = 'black', linewidth = 4)
        plt.xlabel('Time (mins)', fontsize = 20)
        plt.ylabel('Volume', fontsize = 20)
        # plt.title('Alcohol Consumption', fontsize = 15, fontweight = 'bold')
        plt.fill_between(t, 30, where=[(t >= 0) and (t <= 20.1) for t in t], color = 'grey', alpha = 0.75, linewidth = 0.05, label = 'Front-Loading') 
        plt.fill_between(t, 30, where=[(t >= 19.9) and (t <= 120) for t in t], color = 'grey', alpha = 0.25, linewidth = 0.05, label = 'Maintenence') 
        plt.xlim(0,120)
        leg = plt.legend(framealpha = 1, fontsize = 20)
        plt.xticks(fontsize=14) 
        plt.yticks(fontsize=14) 
        plt.show()
     

        if save == True:
            f.savefig('/Users/amyrude/Downloads/graph_'+graph+'_updated.png',transparent=True, dpi = 350)
    if graph == "Param":
        f = plt.figure()
        plt.plot(t, y[7], label="$E_{N}$", linewidth = 3, color = 'red')[0]
        plt.plot(t, abs(y[8]), label='-$Dr_{N}$', linewidth = 3, color = 'blue')[0]
        plt.fill_between(t, 30, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS')
        plt.ylim(0,8)
        # plt.title("DA Modulation of NAc Parameters", fontsize = 20, fontweight = 'bold')
        plt.legend(fontsize = 13)
        plt.xlabel('Time (mins)', fontsize = 13)
        plt.xlim(0,1200)
        plt.show()
        plt.yticks(fontsize = 12)
        plt.xticks(fontsize = 12)

        if save == True:
            f.savefig('/Users/amyrude/Downloads/graph_'+graph+'_updated.png',transparent=True, dpi = 350)

    if graph == 'PFC' or 'Insula' or 'STR' or 'VTA' or 'Cortex' or 'Subcort':
        ax[0].set_xlim(0, 35)
        ax[1].set_xlim(115, 120)
        ax[0].spines['right'].set_visible(False)
        ax[1].spines['left'].set_visible(False)
        ax[1].tick_params(left = False) 
        

        if graph =='STR' or 'Subcort':
            ax[0].set_ylim(0,1.15)
            ax[1].set_ylim(0,1.15)
        if graph == 'Cortex':
            ax[0].set_ylim(-0.1,1)
            ax[1].set_ylim(-0.1,1)
        ax[0].fill_between(t, 30, y2 = -1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        ax[1].fill_between(t, 30, y2 = -1,  where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        d = .015  
        kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
        ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
        ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
        kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
        ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
        ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
        ax[0].set_xlabel('Time (mins)', fontsize = 20, labelpad = 1)
        ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10, fontsize = 20)
        ax[0].tick_params(axis="x", labelsize=14) 
        ax[1].tick_params(axis="x", labelsize=14) 
        ax[0].tick_params(axis="y", labelsize=14) 
        ax[1].tick_params(axis="y", labelsize=14) 
        ax[0].xaxis.set_label_coords(0.6, -0.09)
        ax[1].legend(framealpha = 1, fontsize = 20)
        

        if save == True:
            f.savefig('/Users/amyrude/Downloads/graph_'+graph+'_update.png',transparent=True, dpi = 350)
        

        plt.show()

# ind_plots('Cortex', save=True) #graph: PFC, Insula, STR, VTA, Alc, Param
# ind_plots('Subcort', save=True) #graph: PFC, Insula, STR, VTA, Alc, Param
# ind_plots('Alc', save=True) #graph: PFC, Insula, STR, VTA, Alc, Param


def DA_graphs(F, S, M, L, save = False):
    fail = xppaut_model(t,y0,csTOvta = F, noise=True)['Int']
    low = xppaut_model(t,y0,csTOvta = S, noise=True)['Int']
    medium = xppaut_model(t,y0,csTOvta = M, noise=True)['Int']
    high = xppaut_model(t,y0, csTOvta = L, noise=True)['Int']
    
    fig, ax = plt.subplots(1,2, figsize = (7,5))
    ax[0].plot(t, fail[5],  linewidth = 3, color = 'indigo', alpha = 0.25)
    ax[0].plot(t, low[5],  linewidth = 3,color = 'indigo', alpha = 0.5)
    ax[0].plot(t, medium[5],  linewidth = 3, color = 'indigo', alpha = 0.75)
    ax[0].plot(t, high[5], label = 'VTA', linewidth = 3, color = 'indigo', alpha = 1)
    ax[0].set_xlim(-2,33)
    ax[0].set_ylim(0,1.1)
    # ax[0].set_xlabel('Time (mins)', fontsize =20)
    ax[0].set_ylabel('Firing Rate (Hz)', fontsize =20) 
    ax[0].tick_params(axis="x", labelsize=14) 
    ax[1].tick_params(axis="x", labelsize=14) 
    ax[0].tick_params(axis="y", labelsize=14) 
    ax[1].tick_params(axis="y", labelsize=14)


    ax[1].plot(t, fail[3], linewidth = 3, color = 'chocolate', alpha = 0.25)
    ax[1].plot(t, low[3],  linewidth = 3, color = 'chocolate', alpha = 0.45)
    ax[1].plot(t, medium[3],  linewidth = 3, color = 'chocolate', alpha = 0.75)
    ax[1].plot(t, high[3], label = 'NAc', linewidth = 3, color = 'chocolate', alpha = 1)
    ax[1].set_xlim(-2, 33)
    ax[1].set_ylim(0,1.1)
    ax[1].get_yaxis().set_visible(False)
    ax[0].legend(fontsize = 20)
    ax[1].legend(fontsize = 20, loc = 'upper right')

    plt.subplots_adjust(wspace=0.0, hspace=0.01)
    # fig.supxlabel('Time (mins)', fontsize = 20)
    fig.savefig('/Users/amyrude/Downloads/VTA_NAc.png',transparent=True, dpi = 350)

    # ax[1].set_xlabel('Time (mins)', fontsize ='12')
    # ax[1].set_ylabel('Firing Rate (Hz)', fontsize ='12') 

    fig = plt.figure(figsize = (7,5))
    plt.plot(t, fail[6], linewidth = 3, color = 'black', alpha = 0.25)
    plt.plot(t, low[6], linewidth = 3, color = 'black', alpha = 0.5)
    plt.plot(t, medium[6], linewidth = 3, color = 'black', alpha = 0.75)
    plt.plot(t, high[6], label = 'Alcohol', linewidth = 3, color = 'black', alpha = 1)
    plt.ylim(0, max(high[6])+3)
    plt.xlim(0,120)
    plt.xlabel('Time (mins)', fontsize ='20')
    plt.ylabel('Volume', fontsize ='20')
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.legend(fontsize = 20)
    plt.tight_layout()
    fig.savefig('/Users/amyrude/Downloads/DA_Alc.png',transparent=True, dpi = 350)

    plt.show()

    f, ax = plt.subplots(3,1, figsize = (7,7))
    for k in range(3):
        ax[k].clear()
        ax[k].set_xlim(0,60)
        ax[k].tick_params(axis="x", labelsize=9) 
        ax[k].tick_params(axis="y", labelsize=9) 

    ax[0].get_xaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)

    ax[0].plot(t, fail[5],  linewidth = 3, color = 'indigo', alpha = 0.25)
    ax[0].plot(t, low[5],  linewidth = 3,color = 'indigo', alpha = 0.5)
    ax[0].plot(t, medium[5],  linewidth = 3, color = 'indigo', alpha = 0.75)
    ax[0].plot(t, high[5], label = 'VTA', linewidth = 3, color = 'indigo', alpha = 1)
    ax[1].plot(t, fail[3], linewidth = 3, color = 'chocolate', alpha = 0.25)
    ax[1].plot(t, low[3],  linewidth = 3, color = 'chocolate', alpha = 0.45)
    ax[1].plot(t, medium[3],  linewidth = 3, color = 'chocolate', alpha = 0.75)
    ax[1].plot(t, high[3], label = 'NAc', linewidth = 3, color = 'chocolate', alpha = 1)
    ax[2].plot(t, fail[6], linewidth = 3, color = 'black', alpha = 0.25)
    ax[2].plot(t, low[6],  linewidth = 3, color = 'black', alpha = 0.45)
    ax[2].plot(t, medium[6],  linewidth = 3, color = 'black', alpha = 0.75)
    ax[2].plot(t, high[6], label = 'Alcohol', linewidth = 3, color = 'black', alpha = 1)

    
    ax[0].set_ylim(-0.1,1.1)

    ax[1].set_ylim(-0.1,1.1)
    
    ax[2].set_ylim(-0.1, 27)

    ax[2].set_xlabel('Time (mins)', fontsize = 12)
    ax[0].set_ylabel('Firing Rate (Hz)', fontsize =12) 
    ax[1].set_ylabel('Firing Rate (Hz)', fontsize =12) 
    ax[2].set_ylabel('Volume', fontsize =12) 
    ax[0].legend()        
    ax[1].legend()
    ax[2].legend( loc = 'upper right')

    plt.subplots_adjust(wspace=0, hspace=0.0)
    plt.show()

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

    # f, ax1 = plt.subplots(figsize = (8,6))
    # color = 'tab:red'
    # ax1.set_xlabel('Peak VTA Activity (Hz)', fontsize = '12')
    # ax1.set_ylabel('Total Alcohol Consumed', color=color, fontsize = '12')
    # ax1.scatter(peak_vta, final_alcohol, color=color, linewidths = 0.05)
    # ax1.tick_params(axis='y', labelcolor=color)

    # ax2 = ax1.twinx()  
    # color = 'tab:blue'
    # ax2.set_ylabel('Peak NAc Activity (Hz)', color=color, fontsize = '12', rotation = -90, labelpad = 20)  
    # ax2.scatter(peak_vta, peak_nac, color=color, linewidths = 0.05)
    # ax2.tick_params(axis='y', labelcolor=color)

    # plt.title('Effect of DA Release', fontsize = '15', fontweight = 'bold')
    # fig.tight_layout()  
    # plt.show()

    if save == True:
        fig.savefig('/Users/amyrude/Downloads/DA_FL.png',transparent=True, dpi = 350)
        f.savefig('/Users/amyrude/Downloads/total_plot.png', transparent=True, dpi=350)

# DA_graphs(1.35, 1.5, 1.65, 2, save=True)



def weight_der():
    param_array = np.linspace(1.5,3.5,50)
    der_ratio = []
    f, ax = plt.subplots(1,2)
    for k in range(len(param_array)):
        # init = xppaut_model(t,y0, seekTObin=param_array[k], binTOseek=param_array[k])
        y = der_model(t[208],[0.1, 1, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN], seekTObin=param_array[k], binTOseek=param_array[k] )
        der_ratio.append(y[1]/y[2])
    plt.scatter(param_array, der_ratio, linewidth = 0.005)
    
    plt.show()


def td_vect_ani(t,y0, save = False):
    time = np.array([0,20,40])
    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('Seek', **tfont, fontsize = 15)
    ax.set_ylabel('Binge', **tfont, fontsize = 15)
    ax.set_zlabel('Time',rotation = 90, **tfont, fontsize = 15)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '0.2', '0.4', '0.6', '0.8', '1.0'],  **tfont, fontsize = 12)
    ax.set_zticks([0, 20/200, 40/200], ['0', '20', '40'], rotation=20, **tfont, fontsize = 12)

    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    t = np.linspace(0,time[-1], 500)
    init = xppaut_model(time, y0, csTOvta)['Int'] # Solving the system with IC

    for k in np.arange(len(time)):
        seek_array, binge_array , z = np.meshgrid(np.arange(0, 1.2, 0.2), np.arange(0, 1.2, 0.2), time[k]/200)
        dseek_array = np.zeros(seek_array.shape)
        dbinge_array = np.zeros(binge_array.shape)
        for i in np.arange(len(seek_array)):
            for j in np.arange(len(binge_array)):
                y = [init[0, k], init[1, k], init[2, k], init[3, k], init[4,k], init[5, k], init[6, k], init[7, k], init[8,k]]
                y[1] = seek_array[i, j, 0]
                y[2] = binge_array[i, j, 0]
                deriv = der_model(time[k],y)
                dseek_array[i,j] = deriv[1]
                dbinge_array[i,j] =  deriv[2]
        ax.quiver(seek_array, binge_array, time[k]/200, dseek_array, dbinge_array, 0 ,length = 0.09, alpha = 0.6, arrow_length_ratio = 0.15)
    
    traj_one = xppaut_model(t, y0, csTOvta)['Int']
    ax.scatter3D(traj_one[1][0], traj_one[2][0], t[0]/200, color = 'red', linewidth = 5)
    def update(z):
        ax.plot3D(traj_one[1][:z], traj_one[2][:z], t[:z]/200, color = 'black', linewidth = 1.8)
        return ax
    ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1,repeat=False)
    if save ==True:
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulateNAcparam.gif', writer=writer)
    plt.show()  


def da_animation(param, save=False):
    fig, ax = plt.subplots(3,1, figsize = (7,7))

    def update(z):
        for k in range(3):
            ax[k].clear()
            ax[k].set_xlim(0,60)
            ax[k].tick_params(axis="x", labelsize=9) 
            ax[k].tick_params(axis="y", labelsize=9) 

        ax[0].get_xaxis().set_visible(False)
        ax[1].get_xaxis().set_visible(False)

        y = xppaut_model(t,y0,csTOvta = param[z], noise=False)['Int']        
        ax[0].plot(t, y[5],  linewidth = 3, color = 'indigo', alpha = 1)
        ax[1].plot(t, y[3],  linewidth = 3, color = 'chocolate', alpha = 1)
        ax[2].plot(t, y[6],  linewidth = 3, color = 'black', alpha = 1)
        
        ax[0].set_ylim(-0.1,1.1)

        ax[1].set_ylim(-0.1,1.1)
        
        ax[2].set_ylim(-0.1,27)

        ax[2].set_xlabel('Time (mins)', fontsize = 12)
        ax[0].set_ylabel('Firing Rate (Hz)', fontsize =12) 
        ax[1].set_ylabel('Firing Rate (Hz)', fontsize =12) 
        ax[2].set_ylabel('Volume', fontsize =12) 
        ax[0].legend(['VTA'])        
        ax[1].legend(['NAc'])
        ax[2].legend(['Alcohol'], loc = 'upper right')


        plt.subplots_adjust(wspace=0, hspace=0.0)


        return ax
    ani = animation.FuncAnimation(fig, update, frames=len(param), interval=10,repeat=False)
    if save ==True:
        writer = PillowWriter(fps=30)
        ani.save('/Users/amyrude/Downloads/DAmodulateNAc.gif', writer=writer)
    plt.show()  

   
param_array = np.linspace(1,4,1000)
# da_animation(param_array, save=True)



    
# td_vect_ani(t,y0)    
# td_vect(t,y0)
# ind_plots('Cortex', save=True) #graph: PFC, Insula, STR, VTA, Alc, Param



# #ignore
# import seaborn as sbn
# test = np.array([[1,1,2],[2,2,4]])
# fig = plt.figure()
# # sbn.heatmap(test, cmap = 'Greys')
# plt.plot([0,3,4],[5,6,7])
# plt.xlabel('Low DA    High DA', fontsize = 15)
# plt.show()
# fig.savefig('/Users/amyrude/Downloads/heatmap.png', transparent=False, dpi=350)

def ind_plots_new(csTOvta=csTOvta, t=t, save=False):
    #graph: PFC, Insula, STR, VTA, Alc, Param
    y = xppaut_model(t, y0, noise=True)['Int']
# PFC !!!!!!!!
    f, ax = plt.subplots(1, 2, figsize = (7,5), sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})
    for z in range(2):
        ax[z].plot(t,y[0], label="Setpoint", color = 'royalblue', linewidth = 4)
        ax[z].plot(t, y[1], label="Seek", color = 'midnightblue', linewidth = 4)
        ax[z].plot(t, (y[0]+y[1])/2, '--', label="mPFC Average", color = 'lightblue', linewidth = 4)
        f.suptitle('mPFC Activity', fontsize = 17, fontweight = 'bold')
    ax[0].set_xlim(0, 35)
    ax[1].set_xlim(115, 120)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].tick_params(left = False) 
    
    for k in range(2):
        ax[k].set_ylim(0,1.1)
        ax[k].fill_between(t, 30, y2 = -1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        ax[k].tick_params(axis="x", labelsize=12) 
        ax[k].tick_params(axis="y", labelsize=12) 
    d = .015  
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
    ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
    ax[0].set_xlabel('Time (mins)', fontsize = 15, labelpad = 1)
    ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10, fontsize = 15)
    ax[0].xaxis.set_label_coords(0.6, -0.09)
    ax[1].legend(fontsize = 15) 
    if save == True:
        f.savefig('/Users/amyrude/Downloads/graph_pfc_new.png',transparent=True, dpi = 350)
    
# INSULA!!!!
    f, ax = plt.subplots(1, 2, figsize = (7,5), sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})
    for z in range(2):
        ax[z].plot(t, y[2], label="Binge", color = 'darkgreen', linewidth = 4)

        f.suptitle('Insular Activity', fontsize = 17, fontweight = 'bold')
    ax[0].set_xlim(0, 35)
    ax[1].set_xlim(115, 120)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].tick_params(left = False) 
    
    for k in range(2):
        ax[k].set_ylim(0,1.1)
        ax[k].fill_between(t, 30, y2 = -1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        ax[k].tick_params(axis="x", labelsize=12) 
        ax[k].tick_params(axis="y", labelsize=12) 
    d = .015  
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
    ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
    ax[0].set_xlabel('Time (mins)', fontsize = 15, labelpad = 1)
    ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10, fontsize = 15)
    ax[0].xaxis.set_label_coords(0.6, -0.09)
    ax[1].legend(fontsize = 15) 

    if save == True:
        f.savefig('/Users/amyrude/Downloads/graph_insula_new.png',transparent=True, dpi = 350)
# NAc/ DLS
    f, ax = plt.subplots(1, 2, figsize = (7,5), sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})
    for z in range(2):
        ax[z].plot(t, y[3], label="NAc", color = 'chocolate', linewidth = 4)
        ax[z].plot(t, y[4], label="DLS", color='darkorange', linewidth = 4)
        ax[z].plot(t, y[3]+y[4], '--',label="Striatum", color='tan', linewidth = 4)

        f.suptitle('Striatal Activity', fontsize = 17, fontweight = 'bold')
    ax[0].set_xlim(0, 35)
    ax[1].set_xlim(115, 120)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].tick_params(left = False) 
    
    for k in range(2):
        ax[k].set_ylim(0,1.1)
        ax[k].fill_between(t, 30, y2 = -1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        ax[k].tick_params(axis="x", labelsize=12) 
        ax[k].tick_params(axis="y", labelsize=12) 
    d = .015  
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
    ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
    ax[0].set_xlabel('Time (mins)', fontsize = 15, labelpad = 1)
    ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10, fontsize = 15)
    ax[0].xaxis.set_label_coords(0.6, -0.09)
    ax[1].legend(fontsize = 15) 

    if save == True:
        f.savefig('/Users/amyrude/Downloads/graph_str_new.png',transparent=True, dpi = 350)
    
# VTA !!!!
    f, ax = plt.subplots(1, 2, figsize = (7,5), sharey=True, facecolor='w', gridspec_kw={'width_ratios': [10, 1]})
    for z in range(2):
        ax[z].plot(t, y[5], label="VTA", color = 'indigo', linewidth = 4)
        f.suptitle('VTA Activity', fontsize = 17, fontweight = 'bold')

    ax[0].set_xlim(0, 35)
    ax[1].set_xlim(115, 120)
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].tick_params(left = False) 
    
    for k in range(2):
        ax[k].set_ylim(0,1.1)
        ax[k].fill_between(t, 30, y2 = -1, where=[(t >= 0) and (t <= 3) for t in t], color = 'red', alpha = 0.15, linewidth = 0.05, label = 'CS') 
        ax[k].tick_params(axis="x", labelsize=12) 
        ax[k].tick_params(axis="y", labelsize=12) 
    d = .015  
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d-0.3, +d), (1-d, 1+d), **kwargs)
    ax[1].plot((-d-0.3, +d), (-d, +d), **kwargs)
    ax[0].set_xlabel('Time (mins)', fontsize = 15, labelpad = 1)
    ax[0].set_ylabel('Firing Rate (Hz)', labelpad = 10, fontsize = 15)
    ax[0].xaxis.set_label_coords(0.6, -0.09)
    ax[1].legend(fontsize = 15) 

    if save == True:
        f.savefig('/Users/amyrude/Downloads/graph_vta_new.png',transparent=True, dpi = 350)
# Alcohol !!!!
    f = plt.figure(figsize =(10,7)) #figsize = (14,6)
    plt.plot(t, y[6],  color = 'black', linewidth = 4)
    plt.xlabel('Time (mins)', fontsize = 20)
    plt.ylabel('Volume', fontsize = 20)
    # plt.title('Alcohol Consumption', fontsize = 15, fontweight = 'bold')
    plt.fill_between(t, 30, where=[(t >= 0) and (t <= 20.1) for t in t], color = 'grey', alpha = 0.75, linewidth = 0.05, label = 'Front-Loading') 
    plt.fill_between(t, 30, where=[(t >= 19.9) and (t <= 120) for t in t], color = 'grey', alpha = 0.25, linewidth = 0.05, label = 'Maintenence') 
    plt.xlim(0,120)
    leg = plt.legend(framealpha = 1, fontsize = 20)
    plt.xticks(fontsize=14) 
    plt.yticks(fontsize=14) 

    if save == True:
        f.savefig('/Users/amyrude/Downloads/graph_alc_new.png',transparent=True, dpi = 350)


    plt.show()

    
        
# ind_plots_new(save = True) #graph: PFC, Insula, STR, VTA, Alc, Param


# NULL CLINE ANIMATION

y0 = [0.1, 0, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]
ytraj = [0.1, 0, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]
t= np.linspace(0,50, 1000)
def vector_field(y0, y_traj, t, save=False): 
    #arguments: initial conditions, t, n/m=neuron population numbers (same as IC), name = ['pop one', 'pop two'], save =='yes' if you want to save]
    x1min, x1max, numptsx1, x2min, x2max, numptsx2 = -0.1, 1, 12, -0.1, 1, 12 #setting up the grid space
    x1array, x2array = np.meshgrid(np.linspace(x1min, x1max, numptsx1), np.linspace(x2min, x2max, numptsx2))
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
        for i in np.arange(len(x1array)):
            for j in np.arange(len(x2array)):
                y = [init[0, z], init[1, z], init[2, z], init[3, z], init[4,z], init[5, z], init[6, z], init[7, z], init[8,z]]
                y[1] = x1array[i, j]
                y[2] = x2array[i, j]
                deriv = der_model(t[z],y)
                dx1dt_array[i,j] = deriv[1]
                dx2dt_array[i,j] =  deriv[2]
        ax1.quiver(x1array, x2array, dx1dt_array, dx2dt_array, alpha=0.25, width = 0.0025)    
        
        ax1.plot(init[1, 0:z], init[2, 0:z], color='black', linewidth = 2)
        ax1.plot(init[1,z],init[2,z],marker='o', color = 'red', markersize='10', zorder=10)

        x1list_fine = np.linspace(x1min, x1max, 250)
        x2list_fine = np.linspace(x2min, x2max, 250)
        nullcline = np.zeros((2, 250)) 
        for i in np.arange(250):
            for k in np.arange(9):
                y0[k] = init[k,z]
            y0[1] = x1list_fine[i]
            y0[2] = x2list_fine[i]
        #     #Solving for the nullclines at each time step
            CS = np.heaviside(csDUR-t,0.5)[z] #Conditioned Stimulus
            nullcline[0,i] = F(Eseek * (-spTOseek * y0[0] + csTOseek * CS + binTOseek * y0[2] + seekDRIVE))
            nullcline[1,i] = F(Ebinge * (seekTObin * y0[1] + bingeDRIVE))   
            
        ax1.plot(x1list_fine, nullcline[1, :], 'b-', alpha=0.8, linewidth=1.5)
        ax1.plot(nullcline[0, :], x2list_fine, 'r-', alpha=0.8, linewidth=1.5)
        ax1.set_xlabel('Seek', fontweight='bold', fontsize=15, **tfont)
        ax1.set_ylabel('Binge', fontweight='bold', fontsize=15, **tfont)
        ax1.set_title('Phase Plane Analysis of Seek and Binge', **tfont, fontweight = 'bold', fontsize = '15')

        ax2.plot(t[0:z],init[1,0:z], color = 'red', linewidth = 2)
        ax2.plot(t[0:z],init[2,0:z], color = 'blue', linewidth = 2)
        ax2.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
        ax2.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
        ax2.set_xlim(0,t[-1])
        ax2.set_ylim(0,1.05)        
        # ax2.set_title('Seek and Binge Activity', **tfont, fontweight = 'bold', fontsize = '15')
        ax2.legend(['Seek', 'Binge'])
    
        ax3.plot(t[0:z],init[0,0:z], color = 'darkgreen', label = 'Setpoint', linewidth = 2)
        ax3.set_xlim(0,t[-1])
        ax3.set_ylim(0.09,0.15)
        ax3.set_xlabel('Time', fontweight='bold', fontsize=15, **tfont)
        ax3.set_ylabel('Firing Rate (Hz)', fontweight='bold', fontsize=15, **tfont)
        ax3.set_title('Setpoint Activity', **tfont, fontweight = 'bold', fontsize = '15')
        return ax
    
    ani = animation.FuncAnimation(fig, update, frames=len(t), interval=30,repeat=False)
    if save==True:
        writer = PillowWriter(fps=50)
        ani.save('/Users/amyrude/Downloads/seek_binge_nullani_new.gif', writer=writer)
    plt.show()

# vector_field(y0,ytraj,t, save=True)