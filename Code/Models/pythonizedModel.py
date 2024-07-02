import numpy as np
from params import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
from scipy.integrate import solve_ivp



def derModel(t, y0):
    seek, binge, stop, nac, dls, ALCOHOL = y

    setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
    vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
    aps = Eaps*((1-dlsWeight)*nac + dlsWeight*dls)/2

    dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
    dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
    dstop_dt = (-stop + F(Estop * (binTOstop * binge - spTOstop * setp - stopDRIVE))) / stopTAU
    dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
    ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
    dALCOHOL_dt = (1-dlsWeight)*nac + dlsWeight*dls

    return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

def F(x):
        return 1 / (1 + np.exp(x))
def xppaut_model(t, y0=y0, 
    Ebinge=Ebinge,
    Estop=Estop,
    Enac=Enac,
    Eaps=Eaps,
    Edls=Edls,
    Eseek=Eseek,
    Esetp=Esetp,
    Evta=Evta,

    #TIMESCALES
    seekTAU=seekTAU,
    bingeTAU=bingeTAU,
    stopTAU=stopTAU,
    nacTAU=nacTAU,
    dlsTAU=dlsTAU,

    #DRIVES
    seekDRIVE=seekDRIVE,
    bingeDRIVE=bingeDRIVE,
    stopDRIVE=stopDRIVE,
    nacDRIVE=nacDRIVE,
    dlsDRIVE=dlsDRIVE,

    #SYNAPTIC WEIGHTS
    spTOseek=spTOseek,
    spTOstop=spTOstop,
    seekTOnac=seekTOnac,
    seekTObin=seekTObin,
    binTOnac=binTOnac,
    binTOstop=binTOstop,
    binTOdls=binTOdls,
    stopTObin=stopTObin,
    vtaTOnac=vtaTOnac,
    vtaTOdls=vtaTOdls,
    apsTOseek=apsTOseek,

    #EXTRAS
    dlsWeight=dlsWeight,
    TOLERANCE=TOLERANCE,
    daFACTOR=daFACTOR
):

    def model(t, y):
        seek, binge, stop, nac, dls, ALCOHOL = y

        setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
        aps = Eaps*((1-dlsWeight)*nac + dlsWeight*dls)/2

        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
        dstop_dt = (-stop + F(Estop * (binTOstop * binge - spTOstop * setp - stopDRIVE))) / stopTAU
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
        ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
        dALCOHOL_dt = Eaps*((1-dlsWeight)*nac + dlsWeight*dls)/2

        return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs():
    fig, ax = plt.subplots(figsize=(12, 5))
    y0 = [0, 0.2, 0.2, 0.3, 0, 0]
    t = np.linspace(0, 20, 1000)
    y = xppaut_model(t, y0)
    ALCOHOL = y['Int'][5]
    setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
    vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))

    ax.plot(t, y['Int'][0], label="seek")
    ax.plot(t, y['Int'][1], label="binge")
    ax.plot(t, y['Int'][2], label="stop")
    ax.plot(t, y['Int'][3], label="nac")
    ax.plot(t, y['Int'][4], label="dls")
    ax.plot(t, setp, label="setp")
    ax.plot(t,vta, label="vta")
    ax.plot(t, (y['Int'][3]+y['Int'][4])/2, label="avg")
    #ax.plot(t, ALCOHOL, label="ALCOHOL")
    ax.legend()

    plt.show()


from scipy.optimize import fsolve
def insulaNull(R, binExc = Ebinge, stopExc = Estop):
    binge = [F(binExc * (stopTObin * r  - bingeDRIVE)) / bingeTAU for r in R]
    stop = [F(stopExc * (binTOstop * r  - stopDRIVE)) / stopTAU for r in R]
    return [binge, stop]
def system(inputs, binExc, stopExc):
    return [inputs[0] - F(binExc * (stopTObin * inputs[1]  - bingeDRIVE)) / bingeTAU,
    inputs[1] - F(stopExc * (binTOstop * inputs[0]  - stopDRIVE)) / stopTAU]
def equiFinder(i, j, binExc=Ebinge, stopExc=Estop):
    return fsolve(system, [i, j], args=(binExc, stopExc))
def bingeBif():
    fig, ax = plt.subplots(figsize=(5, 5))
    E = np.linspace(0, 10, 100)
    n=35
    for i in np.linspace(0, 1, n): 
        for j in np.linspace(0, 1, n):
            ax.scatter(E, [equiFinder(i, j, binExc=e)[0] for e in E], c='black', s=.1)
    plt.xlabel("Binge Excitability")
    plt.ylabel("Equilibrium Value")
    plt.savefig("bingeBif", dpi=350)

def stopBif():
    fig, ax = plt.subplots(figsize=(5, 5))
    E = np.linspace(0, 10, 100)
    n=10
    for i in np.linspace(0, 1, n): 
        for j in np.linspace(0, 1, n):
            ax.scatter(E, [equiFinder(i, j, stopExc=e)[1] for e in E], c='black', s=.1)
    plt.xlabel("Stop Excitability")
    plt.ylabel("Equilibrium Value")
    plt.savefig("stopBif", dpi=350)

def binWeightAnim():
    fig, ax = plt.subplots(figsize=(12, 5))
    t = np.linspace(0, 30, 1000)
    y = xppaut_model(t, y0)
    nac = ax.plot(t, y['Int'][3], label="nac")[0]
    dls = ax.plot(t, y['Int'][4], label="dls")[0]
    avg = ax.plot(t, Eaps*(1-dlsWeight)*y['Int'][3]+dlsWeight*y['Int'][4])[0]

    plt.legend()

    def update(frame):
        dlsWeight=frame/frames
        newY = xppaut_model(t, y0, dlsWeight=frame/frames)
        nac.set_data(t, newY['Int'][3])
        dls.set_data(t, newY['Int'][4])
        avg.set_data(t, Eaps*(1-dlsWeight)*newY['Int'][3]+dlsWeight*newY['Int'][4])
        return (nac, dls, avg)

    frames=50
    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=50)
    plt.show()


'''stopBif()
bingeBif()

R = np.linspace(0, 1, 100)
#print([insulaBif(binExc=e) for e in E])
#[insulaBif(binExc=e) for e in np.linspace(0, 10, 100)]
P = {'binExc':2, 'stopExc':Estop}
ax.plot(insulaNull(R, **P)[0], R)
ax.plot(R, insulaNull(R, **P)[1])'''