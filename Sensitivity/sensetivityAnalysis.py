import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *
from SALib.analyze import fast, sobol
from SALib.sample import fast_sampler, saltelli
keys = {
    "Ebinge": 0,
    "Estop": 1,
    "Enac": 2,
    "Eaps": 3,
    "Edls": 4,
    "Eseek": 5,
    # TIMESCALES
    "seekTAU": 6,
    "bingeTAU": 7,
    "stopTAU": 8,
    "nacTAU": 9,
    # DRIVES
    "seekDRIVE": 10,
    "bingeDRIVE": 11,
    "stopDRIVE": 12,
    "nacDRIVE": 13,
    "dlsDRIVE": 14,
    # SYNAPTIC WEIGHTS
    "spTOseek": 15,
    "spSLOPE": 16,
    "spTOstop": 17,
    "seekTOnac": 18,
    "seekTObin": 19,
    "binTOnac": 20,
    "binTOstop": 21,
    "binTOdls": 22,
    "stopTObin": 23,
    "vtaTOnac": 24,
    "vtaTOdls": 25,
    "apsTOseek": 26,
    "TOLERANCE": 27,
    "daFACTOR": 28,
    "vtaSLOPE": 29,
    #INITIAL CONDITIONS
    "seek0":30, 
    "binge0":31, 
    "stop0":32,
    "nac0":33,
    "dls0":34
    }

typeKeys = {
        "S1": 0,
        "ST": 1,
        "S1_conf": 2,
        "ST_conf": 3
    }

def xppaut_model(t,
    Ebinge, 
    Estop, 
    Enac, 
    Eaps, 
    Edls, 
    Eseek,

    #TIMESCALES
    seekTAU, 
    bingeTAU, 
    stopTAU, 
    nacTAU,

    #DRIVES
    seekDRIVE, 
    bingeDRIVE, 
    stopDRIVE, 
    nacDRIVE, 
    dlsDRIVE,

    #SYNAPTIC WEIGHTS
    spTOseek, 
    spTOstop,
    spSLOPE, 
    seekTOnac, 
    seekTObin, 
    binTOnac, 
    binTOstop, 
    binTOdls, 
    stopTObin, 
    vtaTOnac, 
    vtaTOdls, 
    apsTOseek, 

    #EXTRA
    TOLERANCE, 
    daFACTOR, 
    vtaSLOPE,
    
    #INITIAL CONDITIONS
    seek0, binge0, stop0, nac0, dls0
    ):

    y0 = [seek0, binge0, stop0, nac0, dls0, 0]

    def F(x):
        return 1 / (1 + np.exp(x))

    def model(t, y):
        seek, binge, stop, nac, dls, ALCOHOL = y

        setp = 1 - F(spSLOPE * (ALCOHOL - TOLERANCE))
        vta = F(vtaSLOPE*(ALCOHOL - TOLERANCE*daFACTOR))
        aps = Eaps * nac

        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
        dstop_dt = (-stop + F(Estop * (binTOstop * binge  - spTOstop * setp - stopDRIVE))) / stopTAU
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
        ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE)))
        dALCOHOL_dt = (dls + nac) / 2

        return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, 20), y0, dense_output=True)
    y = sol.sol(t)

    return [(y[3, t]+y[4, t])/2 for t in np.arange(1, 1000, 5)]
    return y
    return (y[3, 50]+y[4, 50])/2
    ALCOHOL_integral = np.trapz(sol[:, 3], t)
    return ALCOHOL_integral

def runANALYSIS():
    problem = {
    'num_vars': 35,
    'names': [
    #EXCITABILITY
    'Ebinge','Estop','Enac','Eaps','Edls','Eseek',

    #TIMESCALES'
    'seekTAU','bingeTAU','stopTAU','nacTAU',

    #DRIVES
    'seekDRIVE','bingeDRIVE','stopDRIVE','nacDRIVE','dlsDRIVE',

    #SYNAPTIC WEIGHTS
    'spTOseek','spSLOPE','spTOstop'
    'seekTOnac','seekTObin',
    'binTOnac','binTOstop','binTOdls',
    'stopTObin',
    'vtaTOnac','vtaTOdls',
    'apsTOseek',

    #EXTRA
    'TOLERANCE',
    'daFACTOR',
    'vtaSLOPE',

    #INITIAL CONDITIONS
    'seek0', 'binge0', 'stop0', 'nac0', 'dls0'
    ],

    'bounds': 
    [[0, 10]]*6+
    [[.05, 5]]*4+
    [[-3, 3]]*5+
    [[0, 10]]*12+
    [[.5, 2]]+[[.02, .1]]+[[5, 25]]+
    [[0, 1]]*5
    }

    param_values = fast_sampler.sample(problem, 512)
    t = np.linspace(0, 20, 1000)


    for i in np.arange(0, 200, 1):
        Y = np.array([xppaut_model(t, *params)[i] for params in param_values])
        Si = fast.analyze(problem, Y)
        with open("bigDATAv2.npy", "ab") as f:
            np.save(f, [Si['S1'], Si['ST'], Si['S1_conf'], Si['ST_conf']])

def loadData(names, graph, confidence, file, range=1):
    with open("bigDATAv2.npy", "rb") as f:
        A = np.load(f)

    fig, ax = plt.subplots(figsize=(10, 5))
    n = (np.arange(0, 200, 1))*.1

    for name in names:
        dat = A[:, typeKeys[graph], (keys[name])]
        ax.plot(n, dat, label=str(name))
        if confidence:
            interval = A[:, typeKeys[graph]+2, (keys[name])]
            ax.fill_between(n, dat - interval, dat+interval, alpha=0.3)

    plt.legend()
    plt.ylim(0, range)
    plt.xlim(0, 20)
    plt.xlabel("T")
    if graph=='ST':
        plt.ylabel("Total-Order Sensitivity")
    else:
        plt.ylabel("First-Order Sensitivity")
    #plt.show()
    plt.savefig(str(file)+".png", dpi=250)


def chooseBest(n, graph, lookFor):
    with open("Sensitivity/bigDATAv2.npy", "rb") as f:
        A = np.load(f)
    result = {}
    for name in keys:
        result[name] = simpson(A[:, typeKeys[graph], keys[name]], x=np.linspace(0, 20, 200))
    
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


TAU = ["seekTAU",
    "bingeTAU",
    "stopTAU",
    "nacTAU"]
E = ["Ebinge",
    "Estop",
    "Enac",
    "Eaps",
    "Edls",
    "Eseek"]
DRIVE = ["seekDRIVE",
    "bingeDRIVE",
    "stopDRIVE",
    "nacDRIVE",
    "dlsDRIVE"]
INSULA = ["binTOstop","stopTObin", "Ebinge", "Estop"]
PFC = ["spTOseek", "apsTOseek", "seekTObin", "spTOstop", "seekTOnac"]
INIT = ["seek0","binge0","stop0","nac0","dls0"]
EXTRA = ['TOLERANCE','daFACTOR','vtaSLOPE']
STR_INP = ['vtaTOnac','vtaTOdls', 'binTOnac','binTOdls','seekTOnac']



chooseBest(35, "ST", "DRIVE")
#loadData(["binTOstop", "stopTObin", "seekTObin", "Ebinge", "Estop"], "ST", True, "insula", .3)
#loadData(PFC, "S1", True, "Sensitivity/total/pfc", .35)