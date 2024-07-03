import numpy as np
from params import *
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.animation as animation




def derModel(t, y):
    seek, binge, stop, nac, dls, ALCOHOL = y

    setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
    vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
    aps = Eaps*((1-dlsWeight)*nac + dlsWeight*dls)/2

    dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
    dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
    dstop_dt = (-stop + F(Estop * (binTOstop * binge - spTOstop * setp - stopDRIVE))) / stopTAU
    dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
    ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
    dALCOHOL_dt = ((1-dlsWeight)*nac + dlsWeight*dls)/2

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
        dALCOHOL_dt = ((1-dlsWeight)*nac + dlsWeight*dls)/2

        return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    return {'Int':y, 'Der':[model(t,y0) for t in t]}

def runGraphs():
    fig, ax = plt.subplots(figsize=(12, 5))
    y0 = [0, 0.2, 0.2, 0.3, 0, 0]
    t = np.linspace(0, 50, 1000)
    y = xppaut_model(t, y0, vtaTOdls=0, vtaTOnac=0)
    ALCOHOL = y['Int'][5]
    setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
    vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))

    ax.plot(t, y['Int'][0], label="seek")
    ax.plot(t, y['Int'][1], label="binge")
    ax.plot(t, y['Int'][2], label="stop")
    ax.plot(t, y['Int'][3], label="nac")
    ax.plot(t, y['Int'][4], label="dls")
    ax.plot(t, setp, label="setp")
    #ax.plot(t,vta, label="vta")
    ax.plot(t, (y['Int'][3]+y['Int'][4])/2, label="avg")
    #ax.plot(t, ALCOHOL, label="ALCOHOL")
    ax.legend()

    plt.show()

def insulaNull(R, binExc = Ebinge, stopExc = Estop):
    binge = [F(binExc * (stopTObin * r  - bingeDRIVE)) / bingeTAU for r in R]
    stop = [F(stopExc * (binTOstop * r  - stopDRIVE)) / stopTAU for r in R]
    return [binge, stop]

def insulaBif(binExc = Ebinge, stopExc = Estop):
    A = []
    R = np.linspace(0, 1, 100)
    null = insulaNull(R, binExc, stopExc)
    binge = null[0]
    stop = null[1]
    for i in binge:
         for j in stop:
              if abs(i[0]-j[1])<.007 and abs(i[1]-j[0])<.007:
                   A.append(j)
    
    if len(A)==0:
        return [0, 0, 0]
    A = [A[0], A[len(A)//2], A[-1]]
    print(A)
    return A

def system(inputs, binExc, stopExc):
    return [inputs[0] - F(binExc * (stopTObin * inputs[1]  - bingeDRIVE)) / bingeTAU,
    inputs[1] - F(stopExc * (binTOstop * inputs[0]  - stopDRIVE)) / stopTAU]
def equiFinder(i, j, binExc=Ebinge, stopExc=Estop):
    return fsolve(system, [i, j], args=(binExc, stopExc))

def bingeBif(limit, reso, save=False):
    fig, ax = plt.subplots(figsize=(5, 5))
    E = np.linspace(0, limit, reso)

    for i in np.linspace(.4, .7, limit): 
        for j in np.linspace(.8, 1, limit):
            ax.scatter(E, [equiFinder(i, j, binExc=e)[0] for e in E], c='black', s=.1)
    E = E[math.floor(reso*.295):]
    for i in np.linspace(.4, .7, limit): 
        for j in np.linspace(0, .8, limit):
            ax.scatter(E, [equiFinder(i, j, binExc=e)[0] for e in E], c='black', s=.1)

    plt.xlabel("Binge Excitability")
    plt.ylabel("Equilibrium Value")
    if save:
        plt.savefig("bingeBif", dpi=350)
    else:
        plt.show()

def stopBif(limit, reso, save=False):
    fig, ax = plt.subplots(figsize=(5, 5))
    E = np.linspace(0, limit, reso)
    for i in np.linspace(.5, .6, limit): 
        for j in np.linspace(.4, 1, limit):
            ax.scatter(E, [equiFinder(i, j, stopExc=e)[1] for e in E], c='black', s=.1)
    E = E[math.floor(reso*.44):]
    for i in np.linspace(.5, .6, limit): 
        for j in np.linspace(0, .4, limit):
            ax.scatter(E, [equiFinder(i, j, stopExc=e)[1] for e in E], c='black', s=.1)
    plt.xlabel("Stop Excitability")
    plt.ylabel("Equilibrium Value")
    if save:
        plt.savefig("saveBif", dpi=350)
    else:
        plt.show()



def bothBif(limit, reso, setBin, setStop, save=False):
    fig = plt.figure(figsize=(10, 5))
    t = np.linspace(0, 50, 300)

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    E = np.linspace(0, limit, reso)
    for i in np.linspace(.5, .6, limit): 
        for j in np.linspace(.4, 1, limit):
            ax.scatter([equiFinder(i, j, stopExc=e)[1] for e in E], [equiFinder(i, j, stopExc=e)[0] for e in E], E, c='black', s=.1)
    E = E[math.floor(reso*.44):]
    for i in np.linspace(.5, .6, limit): 
        for j in np.linspace(0, .4, limit):
            ax.scatter([equiFinder(i, j, stopExc=e)[1] for e in E], [equiFinder(i, j, stopExc=e)[0] for e in E], E, c='black', s=.1)
    E = np.linspace(0, limit, reso)
    for ss in setStop:
        y = xppaut_model(t, Estop=ss)['Int']
        binge = y[1]
        stop = y[2]
        ax.scatter(y0[2], y0[1], ss)
        ax.plot(stop, binge, ss, color='red', label='Trajectory')

    ax.view_init(elev=24, azim=63)
    ax.set_xlabel("Stop")
    ax.set_ylabel("Binge")
    ax.set_zlabel("Stop Excitability")

    ax = fig.add_subplot(1, 2, 2, projection='3d')
    for i in np.linspace(.4, .7, limit): 
        for j in np.linspace(.8, 1, limit):
            ax.scatter([equiFinder(i, j, binExc=e)[1] for e in E], [equiFinder(i, j, binExc=e)[0] for e in E], E, c='black', s=.1)
    E = E[math.floor(reso*.295):]
    for i in np.linspace(.4, .7, limit): 
        for j in np.linspace(0, .8, limit):
            ax.scatter([equiFinder(i, j, binExc=e)[1] for e in E], [equiFinder(i, j, binExc=e)[0] for e in E], E, c='black', s=.1)
    ax.plot([],[],[], c='black', label='Equilibrium')
    for sb in setBin:
        y = xppaut_model(t, Ebinge=sb)['Int']
        binge = y[1]
        stop = y[2]
        ax.plot(stop, binge, sb, color='red', label='Trajectory')
        ax.scatter(y0[2], y0[1], sb)
    ax.view_init(elev=25, azim=-140)
    ax.set_xlabel("Stop")
    ax.set_ylabel("Binge")
    ax.set_zlabel("Binge Excitability")

    plt.legend()
    if save:
        plt.savefig("bothBif", dpi=350)    
    else:
        plt.show()

def binWeightAnim(n, save=False):
    fig, ax = plt.subplots(figsize=(12, 5))
    t = np.linspace(0, n, 1000)
    y = xppaut_model(t, y0, vtaTOdls=0, vtaTOnac=0)
    nac = ax.plot(t, y['Int'][3], label="NAc")[0]
    dls = ax.plot(t, y['Int'][4], label="DLS")[0]
    avg = ax.plot(t, Eaps*(1-dlsWeight)*y['Int'][3]+dlsWeight*y['Int'][4], '--', label="Weighted Avg.", color='black')[0]
    def findStop(y):
        for i, s in enumerate(t):
            if (y['Int'][-1][i]) >= TOLERANCE:
                return i, s
            
    flFill = ax.fill_between(t[0:findStop(y)[0]], 0, 1, alpha=0.2, color='purple')
    maFill = ax.fill_between(t[findStop(y)[0]:], 0, 1, alpha=0.2, color='red')
    fl = plt.text(findStop(y)[1]/2.5, .8, 'Front-Loading')
    ma = plt.text(findStop(y)[1]/2.5, .8, 'Maintenance')
    status = plt.text(t[-1]/10, .02, "0% DLS + 100% NAc")

    plt.legend()

    def update(frame):
        dlsWeight=frame/frames
        newY = xppaut_model(t, y0, dlsWeight=frame/frames, vtaTOdls=0, vtaTOnac=0)
        nac.set_data(t, newY['Int'][3])
        dls.set_data(t, newY['Int'][4])
        avg.set_data(t, Eaps*(1-dlsWeight)*newY['Int'][3]+dlsWeight*newY['Int'][4])

        flDummy = ax.fill_between(t[0:findStop(newY)[0]], 0, 1)
        flDp = flDummy.get_paths()[0]
        flDummy.remove()
        flFill.set_paths([flDp.vertices])
        maDummy = ax.fill_between(t[findStop(newY)[0]:], 0, 1)
        maDp = maDummy.get_paths()[0]
        maDummy.remove()
        maFill.set_paths([maDp.vertices])

        fl.set_position((findStop(newY)[1]/2.5, .5))
        ma.set_position(((t[-1]-findStop(newY)[1])/2.5 + findStop(newY)[1], .5))
        status.set_text("Straitum = "+ str(round(dlsWeight*100, 1))+"% DLS + "+str(round(100-dlsWeight*100, 1))+"% NAc")

        return (nac, dls, avg)
    plt.xlabel("T (min)")
    plt.ylabel("Normalized Activity")
    frames=150
    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=10)
    if save:
        FFwriter = animation.FFMpegWriter(fps=60)
        ani.save('changeDLSWeight.mp4', writer = FFwriter)
        plt.close()
    else:
        plt.show()

binWeightAnim(60)
'''stopBif()
bingeBif()

R = np.linspace(0, 1, 100)
#print([insulaBif(binExc=e) for e in E])
#[insulaBif(binExc=e) for e in np.linspace(0, 10, 100)]
P = {'binExc':2, 'stopExc':Estop}
ax.plot(insulaNull(R, **P)[0], R)
ax.plot(R, insulaNull(R, **P)[1])'''


