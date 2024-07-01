import numpy as np
from params import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

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
    return {'Int':y, 'Der':model(t,y0)}

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
    for i, b in enumerate(null[0]):
        for j, s in enumerate(null[1]):
            if abs(i-s) < .01 and abs(j-b) < .01:
                A.append(b)
    if len(A)==0:
        return [1.1, 1.1, 1.1]
    return [A[0], A[len(A)//2], A[-1]]

fig, ax = plt.subplots(figsize=(12, 5))
E = np.linspace(0, 10, 100)
R = np.linspace(0, 1, 100)
#print([insulaBif(binExc=e) for e in E])
ax.plot(E, [insulaBif(binExc=e) for e in E])
'''P = {'binExc':5, 'stopExc':10.5}
ax.plot(insulaNull(R, **P)[0][0], insulaNull(R, **P)[0][1])
ax.plot(insulaNull(R, **P)[1][0], insulaNull(R, **P)[1][1])'''
plt.show()
