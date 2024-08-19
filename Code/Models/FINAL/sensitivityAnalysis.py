import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import *
from SALib.analyze import fast
from SALib.sample import fast_sampler
import os

varNames = {
    'Ebinge': 0, 'Esetp': 1, 'Eseek': 2, 'Evta': 3, 'Edls': 4,
    'seekTAU': 5, 'bingeTAU': 6, 'nacTAU': 7, 'setpTAU': 8, 'vtaTAU': 9, 'dlsTAU': 10,
    'seekDRIVE': 11, 'bingeDRIVE': 12, 'setpDRIVE': 13, 'vtaDRIVE': 14, 'dlsDRIVE': 15,
    'spTOseek': 16, 'seekTOnac': 17, 'seekTObin': 18, 'binTOseek': 19, 'binTOnac': 20, 'vtaTOnac': 21,
    'csTOseek': 22, 'csTOvta': 23, 'csTOdls': 24, 'nacTOsetp': 25, 'dlsTOdls': 26,
    'csDUR': 27, 'EnacDECAY': 28, 'nacDriveDECAY': 29, 'EnacTAU': 30, 'nacDriveTAU': 31,
    'EnacMEAN': 32, 'driveMEAN': 33, 'dlsSCALE': 34
}


typeKeys = {
        "S1": 0,
        "ST": 1,
        "S1_conf": 2,
        "ST_conf": 3
    }


def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

def xppaut_model(t, Ebinge, Esetp, Eseek, Evta, Edls,
seekTAU, bingeTAU, nacTAU, setpTAU, vtaTAU, dlsTAU,
seekDRIVE, bingeDRIVE, setpDRIVE, vtaDRIVE, dlsDRIVE,
spTOseek, seekTOnac, seekTObin, binTOseek, binTOnac, vtaTOnac, csTOseek, csTOvta, csTOdls, nacTOsetp, dlsTOdls,
csDUR, EnacDECAY, nacDriveDECAY, EnacTAU, nacDriveTAU, EnacMEAN, driveMEAN, dlsSCALE
):
    y0 = [0.1, 0, 0.1, 0.1, 0, 0, 0, EnacMEAN, driveMEAN]

    def model(t, y=y0):
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
    return y[3]+y[4]

def runAnalysis():
    os.remove("Code/Models/FINAL/bigDATA.npy")
    problem = {
    'num_vars': 35,
    'names': [
    'Ebinge', 'Esetp', 'Eseek', 'Evta', 'Edls', #5
    'seekTAU', 'bingeTAU', 'nacTAU', 'setpTAU', 'vtaTAU', 'dlsTAU', #6
    'seekDRIVE', 'bingeDRIVE', 'setpDRIVE', 'vtaDRIVE', 'dlsDRIVE', #5
    'spTOseek', 'seekTOnac', 'seekTObin', 'binTOseek', 'binTOnac', 'vtaTOnac', 'csTOseek', 'csTOvta', 'csTOdls', 'nacTOsetp', 'dlsTOdls', #11
    'csDUR', 'EnacDECAY', 'nacDriveDECAY', 'EnacTAU', 'nacDriveTAU', 'EnacMEAN', 'driveMEAN', 'dlsSCALE' #8
    ],

    'bounds': 
    [[.5, 6]]*5+
    [[1, 20]]*6+
    [[-2, 2]]*5+
    [[0, 13]]*11+
    [[0, 3]]+[[0, .01]]*2 + [[.5, 1.5]]*2 + [[.5, 1.5]] + [[-10, 0]] + [[.05, .2]]
    }

    param_values = fast_sampler.sample(problem, 256)
    t = np.linspace(0, 120, reso)

    for i in range(reso):
        Y = np.array([xppaut_model(t, *params)[i] for params in param_values])
        Si = fast.analyze(problem, Y)
        with open("Code/Models/FINAL/bigDATA.npy", "ab") as f:
            np.save(f, [Si['S1'], Si['ST'], Si['S1_conf'], Si['ST_conf']])
        print(str((i/reso)*100)+"%")

def recover():
    arrays = []
    with open("Code/Models/FINAL/bigDATA.npy", "rb") as f:        
        while True:
            try:
                arrays.append(np.load(f))
            except EOFError:
                break
    with open("Code/Models/FINAL/bigDATA.npy", "wb") as f:
        np.save(f, arrays)
    return arrays

def loadData(names, graph, confidence=False, range=1):
    with open("bigDATA.npy", "rb") as f:
        A = np.load(f)

    fig, ax = plt.subplots(figsize=(10, 5))
    t = np.linspace(0, 120, reso)

    for name in names:
        dat = A[:, typeKeys[graph], (varNames[name])]
        ax.plot(t, dat, label=str(name))
        if confidence:
            interval = A[:, typeKeys[graph]+2, (varNames[name])]
            ax.fill_between(t, dat - interval, dat+interval, alpha=0.3)

    plt.legend()
    plt.ylim(0, range)
    plt.xlabel("T")
    if graph=='ST':
        plt.ylabel("Total-Order Sensitivity")
    else:
        plt.ylabel("First-Order Sensitivity")
    
    plt.show()

def chooseBest(n, graph, lookFor):
    with open("Code/Models/FINAL/bigDATA.npy", "rb") as f:
        A = np.load(f)
    result = {}
    for name in varNames:
        result[name] = simpson(A[:, typeKeys[graph], varNames[name]], x=np.linspace(0, 20, 200))
    
    names = list(result.keys())
    values = list(result.values())
    indexes = np.argsort(values)[-n:]
    result = {names[i]: values[i] for i in indexes}

    i=0
    for r in result:
        if lookFor in r:
            print(r+": "+str(result[r]))
            i +=1
            if i==n:
                return

reso = 200

E = ['Ebinge', 'Esetp', 'Eseek', 'Evta', 'Edls']
TAU = ['seekTAU', 'bingeTAU', 'nacTAU', 'setpTAU', 'vtaTAU', 'dlsTAU']
SYN = ['spTOseek', 'seekTOnac', 'seekTObin', 'binTOseek', 'binTOnac', 'vtaTOnac', 'csTOseek', 'csTOvta', 'csTOdls', 'nacTOsetp', 'dlsTOdls']
MISC = ['csDUR', 'EnacDECAY', 'nacDriveDECAY', 'EnacTAU', 'nacDriveTAU', 'EnacMEAN', 'driveMEAN', 'dlsSCALE']


#TO SEE GRAPHS
#loadData(ArrayLike parameterToSee, String 'S1' or 'ST', Boolean showConfidenceIntervals=False, Float yAxisLimit)
#e.g. 
loadData(['seekTAU', 'bingeTAU'], 'S1', confidence=False, range=.1)

#TO SEARCH FOR MOST/LEAST SENSITIVE PARAMETERS
#chooseBest(int numberOfParams, String 'S1' or 'ST', String whatToSearch)
#e.g. chooseBest(30, 'ST', 'vta')