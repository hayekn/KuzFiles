#INIT
y0 = [0, .2, .2, 0.3, 0, 0, 0]

#EXCITABILITIES

Ebinge = 2
Enac = 1
Eaps = 1
Eseek=1
Esetp = 3
Evta = 1

#TIMESCALES
seekTAU = 1
bingeTAU = 1
nacTAU = 1
setpTAU = 1

#DRIVES
seekDRIVE = 0
bingeDRIVE = -2
nacDRIVE = -1.4
vtaDRIVE = -2

#SYNAPTIC WEIGHTS
spTOseek = 10
seekTOnac = 2
seekTObin = 3
binTOnac = 3
binTObin = 1
vtaTOnac = 1
apsTOseek = 1
csTOvta = 3
csTOseek = 2

#NEGATIVE STIM
Ens = 10
nsLEVEL = 1
nsSTART=0
nsDURATION=10
nsTOvta = 3
nsTOseek = 3
nsTObin = 0

#EXTRAS
TOLERANCE = 20
daFACTOR = 0.2
decayFac = .001