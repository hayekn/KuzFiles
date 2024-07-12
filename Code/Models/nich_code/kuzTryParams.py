#setp, seek, binge, nac, vta, av, ALCOHOL

y0 = [0, .5, .1, 0.2, 0, 0, 0]

#EXCITABILITIES

Ebinge = 1
Enac = 1.84
#eav = 1.84
Eav = 1
Eseek = 2
Esetp = 5
Evta = 12
Edls = 2

#TIMESCALES
seekTAU = 1
bingeTAU = .5
nacTAU = 1
setpTAU = 1
vtaTAU = 1
avTAU = 20
dlsTAU = 1

#DRIVES
seekDRIVE = 0
bingeDRIVE = -1.4
nacDRIVE = -1.2
vtaDRIVE = -1.4
setpDRIVE = 0
avDRIVE = -1.4
dlsDRIVE = -5

#SYNAPTIC WEIGHTS
#spTOseek = 5
spTOseek = 5
seekTOnac = 10
seekTObin = 3
seekTOseek = 4.2
seekTOdls = 3
binTOnac = 1
binTObin = 0
avTOseek = 0
avTOsetp = .2
csTOvta = 3
#cstoseek = 6
csTOseek = 3
nacTOav = 1
vtaTOnac = 1
dlsTOdls = 1

#NEGATIVE STIM
Ens = 10
nsLEVEL = 1
nsSTART=0
nsDURATION=10
nsTOvta = 3
nsTOseek = 3
nsTOnseek = 1
nsTObin = 0

#EXTRAS
TOLERANCE = 20
daFACTOR = 0.2
decayFac = .001
setpDELTA = 1
csDUR = 3