import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

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
csTOvta = 2.5 # modulate this connection to change magnitude of DA peak (3 to 5.5)
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

def xppaut_model(t, fuzz=False, csTOvta=csTOvta, y=y0):
    def model(t, y):
        setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE = y
        if fuzz:
            noise = np.random.normal(0,.3)
        else:
            noise = 0

        cs = np.heaviside(csDUR-t, .5)

        dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN-Enac)
        dnacDRIVE_dt =  vta/nacDriveTAU + nacDriveDECAY*(driveMEAN-nacDRIVE)
        
        dsetp_dt = (-setp + F(Esetp * (nacTOsetp * (nac+dlsSCALE*dls) + setpDRIVE)) + noise) / setpTAU
        dseek_dt = (-seek + F(Eseek * (-spTOseek * setp + csTOseek * cs + binTOseek * binge + seekDRIVE)) + noise) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek + bingeDRIVE)) + noise) / bingeTAU
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE)) + noise) / nacTAU
        ddls_dt = (-dls + F(Edls * (dlsTOdls * dls + csTOdls * cs + dlsDRIVE + noise))) / dlsTAU
        
        dvta_dt = (-vta + F(Evta*(csTOvta * cs + vtaDRIVE)) + noise) / vtaTAU
        dALCOHOL_dt = (nac+dlsSCALE*dls)
        return [dsetp_dt, dseek_dt, dbinge_dt, dnac_dt, ddls_dt, dvta_dt, dALCOHOL_dt, dEnac_dt, dnacDRIVE_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    y[4] = dlsSCALE*y[4]
    global cs
    cs = np.heaviside(csDUR-t, .5)
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'CS':cs}


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
    plt.show()


# runGraphs(120)

