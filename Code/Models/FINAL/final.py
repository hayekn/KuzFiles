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
csTOvta = 3 # modulate this connection to change magnitude of DA peak (3 to 5.5)
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
t = np.linspace(0,120,5)
#setp, seek, binge, nac, dls, vta, ALCOHOL, Enac, nacDRIVE

def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

def xppaut_model(t, fuzz, csTOvta=csTOvta):
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



