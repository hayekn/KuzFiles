#INIT
y0 = [0, .2, .2, 0.3, 0, 0, 0]

#EXCITABILITIES
Ebinge = 1
Estop = 7
Enac = 1.84
Eaps = 1
Edls = 1.84
Eseek=1
Esetp = 3
Evta = 6

#TIMESCALES
seekTAU = 1
bingeTAU = 1
stopTAU = 1
nacTAU = 1
dlsTAU = .005

#DRIVES
seekDRIVE = 0
bingeDRIVE = -2
stopDRIVE = 0.5
nacDRIVE = -1.4
dlsDRIVE = -.8

#SYNAPTIC WEIGHTS
spTOseek = 10
spTOstop = 1
spDURATION = 1
spTObin = 3
seekTOnac = 2
seekTObin = 1
seekTOseek = 1
binTOnac = 3
binTOstop = 1
binTOdls = 2
binTObin=1
stopTObin = 2.5
vtaTOnac = 1
vtaTOdls = 1
vtaTObinge = 3
apsTOseek = 1

#NEGATIVE STIM
Ens = 15
nsLEVEL = 1
nsSTART=0
nsDURATION=10
nsTOvta = 5
nsTOstop = 2
nsTOseek = 1
nsTObinge = 5

#EXTRAS
dlsWeight = .3
TOLERANCE = 20
daFACTOR = 0.2
decayFac = .001
csTOvta = 1
csTOseek = 2

#import pythonizedModel
#pythonizedModel