import numpy as np
import params
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from SALib.analyze import fast, sobol
from SALib.sample import fast_sampler, saltelli

def xppaut_model(t, **params.paramDict()):

    #INITIAL CONDITIONS
    y0 = [0, 0.2, 0.2, 0.3, 0, 0]

    def F(x):
        return 1 / (1 + np.exp(x))

    def model(t, y):
        seek, binge, stop, nac, dls, ALCOHOL = y

        setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
        vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))
        aps = Eaps*(nac + dls)/2

        dseek_dt = (-seek + F(Eseek * (spTOseek * setp - apsTOseek * aps - seekDRIVE))) / seekTAU
        dbinge_dt = (-binge + F(Ebinge * (stopTObin * stop - seekTObin * seek - bingeDRIVE))) / bingeTAU
        dstop_dt = (-stop + F(Estop * (binTOstop * binge - spTOstop * setp - stopDRIVE))) / stopTAU
        dnac_dt = (-nac + F(Enac * (-vtaTOnac * vta - seekTOnac * seek - binTOnac * binge - nacDRIVE))) / nacTAU
        ddls_dt = (-dls + F(Edls * (-binTOdls * binge - vtaTOdls * vta - dlsDRIVE))) / dlsTAU
        dALCOHOL_dt = (dls + nac) / 2

        return [dseek_dt, dbinge_dt, dstop_dt, dnac_dt, ddls_dt, dALCOHOL_dt]

    sol = solve_ivp(model, (0, 20), y0, dense_output=True)
    y = sol.sol(t)
    return y


fig, ax = plt.subplots(figsize=(12, 5))
y = xppaut_model(np.linspace(0, 20, 1000))
ALCOHOL = y[5]
setp = 1 - F(Esetp * (ALCOHOL - TOLERANCE))
vta = F(Evta*(ALCOHOL - TOLERANCE*daFACTOR))

'''ax.plot(t, y[0], label="seek")
ax.plot(t, y[1], label="binge")
ax.plot(t, y[2], label="stop")
ax.plot(t, y[3], label="nac")
ax.plot(t, y[4], label="dls")
ax.plot(t, setp, label="setp")
ax.plot(t,vta, label="vta")
ax.plot(t, (y[3]+y[4])/2, label="avg")
ax.plot(t, ALCOHOL, label="ALCOHOL")
ax.legend()'''

plt.show()