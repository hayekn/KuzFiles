#seek:0, nseek:1, seekD1:2, seekD2:3, nseekD1:4, nseekD2:5, 
#seekGPE:6, seekGPI:7, seekSTN:8, nseekGPE:9, nseekGPI:10, nseekSTN:11,  
#binge:12, nac:13, ALCOHOL:14, aps:15, setp:16, vta:17

#INIT
y0 = [.2, .2, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0,
      .2, .2, 0, 0, 0, 0]

#EXCITABILITIES

Ebinge = 2
Enac = 1
Eaps = 1
Eseek = 6
Enseek = 6
Esetp = 3
Evta = 1

#TIMESCALES
seekTAU = 1
nseekTAU = 1
bingeTAU = 1
nacTAU = 1
setpTAU = 1

#DRIVES
seekDRIVE = .6
nseekDRIVE = .6
bingeDRIVE = -2
nacDRIVE = -1.4
vtaDRIVE = -2

#SYNAPTIC WEIGHTS
spTOseek = 1
spTOnseek = 1
seekTOnac = 2
seekTObin = 3
seekTOnseek = 1
nseekTOseek = 1
binTOnac = 3
binTObin = 1
vtaTOnac = 1
apsTOseek = 0.1
csTOvta = 3
csTOseek = .1

#NEGATIVE STIM
Ens = 10
nsLEVEL = 1
nsSTART=0
nsDURATION=10
nsTOvta = 3
nsTOseek = 3
nsTOnseek = 1
nsTObin = 0

#NEW STUFF
seekTOseekD1=1
seekTOseekD2=1
nseekTOnseekD1=1
nseekTOnseekD2=1

#EXTRAS
TOLERANCE = 20
daFACTOR = 0.2
decayFac = .001
setpDELTA = 1