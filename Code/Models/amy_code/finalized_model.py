import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
import random
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
tfont = {'fontname':'Times New Roman'}

#Specify IC and time step

y0 = [0.2, 0, 0.1, 0.1, 0, 0, 0, 0.6, -7] #seek, setp, binge,  nac, dls, ALCOHOL, vta, Enac, Nac Drive
y_traj = [0.2, 0, 0.1, 0.1, 0, 0, 0, 0.6, -7] #seek, setp, binge, nac, dls, ALCOHOL, vta, Enac, Nac Drive; used in phase plane analysis
t= np.linspace(0,30,500)  

## Defining the Parameters
#EXCITABILITIES
Ebinge = 5
Esetp = 1.5
Eseek= 5
Evta = 2
Edls = 3

#TIMESCALES
seekTAU = 1
bingeTAU = 1
nacTAU = 1
setpTAU = 30
vtaTAU = 1
dlsTAU = 1

#DRIVES
seekDRIVE = 1
bingeDRIVE = 1
setpDRIVE = 1
vtaDRIVE = 3.5
dlsDRIVE = 1.4

#SYNAPTIC WEIGHTS
spTOseek = 8
seekTOnac = 2.5
seekTObin = 1.75
binTOseek = 1.75
binTOnac = 2.5
vtaTOnac = 2.5
csTOseek = 5
csTOvta = 4.2 # modulate this connection to change magnitude of DA peak (3 to 5.5)
csTOdls = 2


#EXTRAS
csDUR = 3
decayFac = 0.001
nacWEIGHT = 0.75

#DA Modulation
EnacDECAY = 0.003
nacdriveDECAY = 0.003
EnacTAU = 0.8
nacdriveTAU = 0.8
EnacMEAN = 0.6
driveMEAN = -7
nacdrSCALE = 5
dlsSCALE = 0.1

#Sigmoidal Function
def F(x): # + = excitatory, - = inhibitory
        return 1 / (1 + np.exp(-x))

#Main Model
def binge_model(t, y0, param):
    def model(t, y):
        seek, setp, binge, nac, dls,  ALCOHOL, vta, Enac, nacDRIVE = y

        csTOvta = param # for animation when varying DA concentration

        dEnac_dt = vta/EnacTAU + EnacDECAY*(EnacMEAN - Enac)
        dnacDRIVE_dt = vta/nacdriveTAU + nacdriveDECAY*(driveMEAN - nacDRIVE)

        CS = np.heaviside(csDUR-t, 0.5) #Conditioned Stimulus
        dseek_dt = (-seek + F(Eseek * (binTOseek * binge + csTOseek * CS - spTOseek * setp - seekDRIVE))) / seekTAU #Seek Activity
        dsetp_dt = (-setp + F(Esetp * (0.9*nac + dlsSCALE * dls - setpDRIVE - setpDRIVE))) / setpTAU #Alcohol Variable
        dbinge_dt = (-binge + F(Ebinge * (seekTObin * seek - bingeDRIVE))) / bingeTAU #Binge Activity
        dnac_dt = (-nac + F(Enac * (vtaTOnac * vta + seekTOnac * seek + binTOnac * binge + nacDRIVE))) / nacTAU #NAc Activity
        ddls_dt = (-dls + F(Edls * ( 3*dls + csTOdls * CS - dlsDRIVE)))/ dlsTAU
        dALCOHOL_dt = nac + dlsSCALE * dls # Alcohol consumed 
        dvta_dt = (-vta + F(Evta*( csTOvta * CS - vtaDRIVE))) / vtaTAU #VTA activity
       

        return [dseek_dt, dsetp_dt, dbinge_dt, dnac_dt, ddls_dt, dALCOHOL_dt, dvta_dt, dEnac_dt, dnacDRIVE_dt]

    sol = solve_ivp(model, (0, t[-1]), y0, dense_output=True)
    y = sol.sol(t)
    CS = np.heaviside(csDUR-t,0.5)
    y[4] = dlsSCALE * y[4]
    return {'Int':y, 'Der':[model(t,y0) for t in t], 'CS':CS}
