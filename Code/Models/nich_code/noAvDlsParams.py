#setp, seek, binge, nac, vta, av, ALCOHOL

#EXCITABILITIES

Ebinge = 2
#Enac = .8
#eav = 1.84
Eseek = 2
Esetp = 1.84
Evta = 12
Edls = 2

#TIMESCALES
seekTAU = 1
bingeTAU = .5
nacTAU = 1
setpTAU = 25
vtaTAU = 1
dlsTAU = 1

#DRIVES
seekDRIVE = 0
bingeDRIVE = -1
nacDRIVE = -4
vtaDRIVE = -1.4
setpDRIVE = -1
dlsDRIVE = -1.5

#SYNAPTIC WEIGHTS
#spTOseek = 5
spTOseek = 10
seekTOnac = 4
seekTObin = 2.5
seekTOseek = 4
seekTOdls = 3
binTOnac = 1
binTOseek = 2.5
csTOvta = 3
#cstoseek = 6
csTOdls = 1
csTOseek = 4
nacTOsetp = .5 #BINGE TIME
vtaTOnac = 0
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

#PROX
driveMEAN = -.5
proxMEAN = .2
proxDECAY = .005
proxTAU = 2