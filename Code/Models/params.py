Ebinge = 10.5
Estop = 10.5
Enac = 1.84
Eaps = 1
Edls = 1.84
Eseek=1
Esetp = 6
Evta = 12

#TIMESCALES
seekTAU = 1
bingeTAU = 1
stopTAU = 1
nacTAU = 1
dlsTAU = 1

#DRIVES
seekDRIVE = 0.01
bingeDRIVE = 0.5
stopDRIVE = 0.5
nacDRIVE = -1.4
dlsDRIVE = -1.4

#SYNAPTIC WEIGHTS
spTOseek = 5
spTOstop = 1
seekTOnac = 10
seekTObin = 3
binTOnac = 1
binTOstop = 1
binTOdls = 2.5
stopTObin = 2
vtaTOnac = 1
vtaTOdls = 1
apsTOseek = 1

#EXTRAS
dlsWeight = 0.5
TOLERANCE = 5
daFACTOR = 0.1


def paramDict():
    D = {

    #EXCITABILITIES  
    "Ebinge": 10.5,
    "Estop": 10.5,
    "Enac": 1.84,
    "Eaps": 1,
    "Edls": 1.84,
    "Eseek": 1,
    "Esetp": 6,
    "Evta": 12,

    #TIMESCALES
    "seekTAU": 1,
    "bingeTAU": 1,
    "stopTAU": 1,
    "nacTAU": 1,
    "dlsTAU": 1,

    #DRIVES
    "seekDRIVE": 0.01,
    "bingeDRIVE": 0.5,
    "stopDRIVE": 0.5,
    "nacDRIVE": -1.4,
    "dlsDRIVE": -1.4,

    #SYNAPTIC WEIGHTS
    "spTOseek": 5,
    "spTOstop": 1,
    "seekTOnac": 10,
    "seekTObin": 3,
    "binTOnac": 1,
    "binTOstop": 1,
    "binTOdls": 2.5,
    "stopTObin": 2,
    "vtaTOnac": 1,
    "vtaTOdls": 1,
    "apsTOseek": 1,

    #EXTRAS
    "TOLERANCE": 10,
    "daFACTOR": 0.1
    
    }

    return D