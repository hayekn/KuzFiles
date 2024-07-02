import numpy as np
from params import *
import matplotlib.pyplot as plt
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
    dALCOHOL_dt = ((1-dlsWeight)*nac + dlsWeight*dls)/2

    return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

def F(x):
        return 1 / (1 + np.exp(x))
def xppaut_model(t, y0):

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

from scipy.optimize import fsolve
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


'''stopBif()
bingeBif()

R = np.linspace(0, 1, 100)
#print([insulaBif(binExc=e) for e in E])
#[insulaBif(binExc=e) for e in np.linspace(0, 10, 100)]
P = {'binExc':2, 'stopExc':Estop}
ax.plot(insulaNull(R, **P)[0], R)
ax.plot(R, insulaNull(R, **P)[1])'''


print(xppaut_model(np.linspace(0, 20, 100), [0, 0, 0, 0, 0, 0])['Der'])