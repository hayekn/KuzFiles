import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast, sobol
from SALib.sample import fast_sampler, saltelli

# Define the model function
def xppaut_model(t,
    #EXCITABILITIES
    Ebinge = 10.5,
    Estop = 10.5,
    Enac = 1,
    Eaps = 1,
    Edls = 1,
    Eseek=1,

    #TIMESCALES
    seekTAU = 1,
    bingeTAU = 1,
    stopTAU = 1,
    nacTAU = 1,

    #DRIVES
    seekDRIVE = 0.01,
    bingeDRIVE = 0.01,
    stopDRIVE = 0.01,
    nacDRIVE = -1,
    dlsDRIVE = 0.01,

    #SYNAPTIC WEIGHTS
    spTOseek = 5,
    spSLOPE = 6,
    seekTOnac = 10,
    seekTObin = 3,
    binTOnac = 1,
    binTOstop = 1,
    binTOdls = 2.5,
    stopTObin = 2,
    vtaTOnac = 1,
    vtaTOdls = 1,
    apsTOseek = 1,
    nsSTART = 0,
    nsDUR = 0,
    nsTOstop = 0,
    nacCAP = 1,
    TOLERANCE = 5,
    daFACTOR = 0.25,
    vtaSLOPE = 4):

    #INITIAL CONDITIONS
    y0 = [0, 0.2, 0.2, 0.3, 0, 0]

    def F(x):
        return 1 / (1 + np.exp(x))

    def model(t, y):
        seek, binge, stop, nac, dls, ALCOHOL = y

        setp = 1 - F(spSLOPE * (ALCOHOL - TOLERANCE))
        vta = F(vtaSLOPE*(ALCOHOL - TOLERANCE*daFACTOR))
        aps = Eaps * nac

        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
        dstop_dt = (-stop + F(Estop * (binTOstop * binge  - stopDRIVE))) / stopTAU
        dnac_dt = (-nac / nacCAP + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
        ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE)))
        dALCOHOL_dt = (dls + nac) / 2

        return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

    # Solve ODEs
    sol = solve_ivp(model, (0, 20), y0, dense_output=True)
    y = sol.sol(t)

    '''fig, ax = plt.subplots(figsize=(12, 5))
    ALCOHOL = y[5]
    setp = 1 - F(spSLOPE * (ALCOHOL - TOLERANCE))
    ax.plot(t, setp)
    plt.show()'''

    return [y[-1, t] for t in np.arange(1, 1000, 10)]
    return [(y[3, t]+y[4, t])/2 for t in np.arange(1, 1000, 10)]
    return (y[3, 50]+y[4, 50])/2
    ALCOHOL_integral = np.trapz(sol[:, 3], t)
    return ALCOHOL_integral


'''problem = {
    'num_vars': 1,
    'names': ['Test'],
    'bounds': [[.1, 5]]
    }
param_values = fast_sampler.sample(problem, 128)
t = np.linspace(0, 20, 100)
Y = np.array([xppaut_model(t, vtaTOdls=sample[0]) for sample in param_values])
Si = fast.analyze(problem, Y)
print(Si['S1'])'''


def runANALYSIS(steps, range, name):
    problem = {
    'num_vars': 1,
    'names': ['NA'],
    'bounds': [range]
    }

    param_values = fast_sampler.sample(problem, 128)
    t = np.linspace(0, 20, 1000)

    match name:
        case "spTOseek":
            Z = np.array([[xppaut_model(t, spTOseek=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "apsTOsp":
            Z = np.array([[xppaut_model(t, apsTOsp=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Ebinge":
            Z = np.array([[xppaut_model(t, Ebinge=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Estop":
            Z = np.array([[xppaut_model(t, Estop=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Eseek":
            Z = np.array([[xppaut_model(t, Eseek=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Eaps":
            Z = np.array([[xppaut_model(t, Eaps=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Enac":
            Z = np.array([[xppaut_model(t, Enac=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "Edls":
            Z = np.array([[xppaut_model(t, Edls=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "vtaTOnac":
            Z = np.array([[xppaut_model(t, vtaTOnac=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "nacTAU":
            Z = np.array([[xppaut_model(t, nacTAU=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "seekTAU":
            Z = np.array([[xppaut_model(t, seekTAU=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "stopTAU":
            Z = np.array([[xppaut_model(t, stopTAU=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])
        case "bingeTAU":
            Z = np.array([[xppaut_model(t, bingeTAU=sample[0])[i] for sample in param_values] for i in np.arange(0, 100, steps)])

   
    A = []
    for z in Z:
        resultA = z
        Si = fast.analyze(problem, resultA)
        A.append(Si['S1'][0])

    n = (np.arange(0, 100, steps))*.2
    ax.plot(n, A, label=str(name))


    plt.ylim(0, 1)
    plt.xlim(0, 20)
    plt.xlabel("T")
    plt.ylabel("Normalized Sensitivity")


fig, ax = plt.subplots(figsize=(12, 5))
E = ["spTOseek"]
for e in E:
    runANALYSIS(2, [0, 12], str(e))
plt.legend()
plt.show()



'''
Y = np.array([xppaut_model(sample[0], t)[-1] for sample in param_values])

# Perform the sensitivity analysis
Si = fast.analyze(problem, Y)

# Print the results
print("First-order sensitivity indices:", Si['S1'])
print("Total-order sensitivity indices:", Si['ST'])'''
