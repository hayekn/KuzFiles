#INIT
y0 = [0, .2, .2, 0.3, 0, 0]

#EXCITABILITIES
Ebinge = 7
Estop = 7
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
seekTOnac = 4
seekTObin = 3
binTOnac = 1
binTOstop = 1
binTOdls = 1
stopTObin = 2.5
vtaTOnac = 1
vtaTOdls = 1
apsTOseek = 1

#NEGATIVE STIM
Ens = 2
nsLEVEL = 1
nsSTART=15
nsDURATION=10
nsTOvta = 5
nsTOstop = 1

#EXTRAS
dlsWeight = .3
TOLERANCE = 20
daFACTOR = 0.1

import pythonizedModel
pythonizedModel