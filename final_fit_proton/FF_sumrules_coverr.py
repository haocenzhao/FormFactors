######## 
### FF_sumrules_coverr.py
### Produces errors on form factors with sum rules at fixed Q^2 points.
### Accepts five arguments from command line:
### 1. Q2max points for fitting data, separated by commas.
### 2 .dat filename of completed fit.
### 3. Filename for output.
########

import sys # for user options from command line
from math import *
from numpy import array, dot
import numpy as np

########################
### CONSTANTS
########################
pi = 3.14159265358
alpha = 7.29735e-3 # Fine-structure constant.
me = 5.10999e-4 # Mass of the electron in GeV.
Mp = 0.938272 # Mass of the proton in GeV.
mpi = 0.13957 # Mass of the charged pion in GeV.
GE0 = 1. # G_E(0), value of the Sachs electric form factor at q^2 = 0, which for the proton is 1.
GM0 = 2.792847356 # G_M(0), value of the Sachs magnetic form factor at q^2 = 0, which for the proton is the magnetic moment.
Lambda2 = 0.71 # Best fit value for parameter in dipole functional form.
GeVfm = 0.197327 # GeV^(-1) to fm conversion: 1 GeV^(-1) = 0.197327 fm.
barn = 0.389379338e-3 # 1 GeV^(-2)=0.389e-3 barn.
########################

# Input.
Q2max = float(sys.argv[1])
fitdata = open(sys.argv[2], 'r').readlines()
outfile = sys.argv[3]
fittype = sys.argv[4]

logQ2min = -6.0 #1e-6
logQ2max = 4.0  #1e4
logQ2step = 0.005
#logQ2step = 0.01
logN = int((logQ2max-logQ2min)/logQ2step)
Q2 = np.zeros(logN, dtype=float)
for i in range(logN):
   logQ2 = logQ2min + i*logQ2step
   Q2[i] = 10.0**logQ2

# Q2int = 0.3
# Q2step = 0.0005
# Q2 = [Q2step*i for i in range(0,int(Q2int/Q2step)+1)]
# Q2fix=100
# Q2step = 0.05
# Q2 = Q2+[Q2int+Q2step*i for i in range(int((Q2fix-Q2int)/Q2step)+2)]
# Q2step = 1.
# Q2 = Q2+[Q2fix+Q2step*i for i in range(int((Q2max-Q2fix)/Q2step)+2)]


# Extract relevant fit data.
[Q2max, tcut, t0] = [float(fitdata[0].split()[i]) for i in range(0,3)]
[ndof, num_tot, num_Mainz, num_world, num_pol, num_fake, kmax_temp, kmax_temp, kslope, kMainz, kworld, num_fakeHQ, num_SBS] = [float(fitdata[1].split()[i]) for i in range(0,13)]
[redchi2, chi2, chi2xsMainz, chi2xsworld, chi2pol, chi2gefake, chi2gmfake, chi2gecoef, chi2gmcoef, chi2csyst, chi2NMainz, chi2Nworld, chi2gefakeHQ, chi2gmfakeHQ, chi2SBS, chi2erad, chi2mrad
] = [float(fitdata[2].split()[i]) for i in range(0,17)]

[erad, derad, mrad, dmrad] = [float(fitdata[3].split()[i]) for i in range(0,4)]

gecoef_fit = [float(x) for x in fitdata[4].lstrip('[').rstrip(']\n').split(',')]
gmcoef_fit = [float(x) for x in fitdata[5].lstrip('[').rstrip(']\n').split(',')]

kmax = len(gecoef_fit)-1
nmax = kmax-4

if fittype==3:
    NMainz = np.array([float(x) for x in fitdata[7].lstrip('[').rstrip(']\n').split(',')])
else:
    NMainz = np.zeros(1, dtype=float)

if fittype==2:
    Nworld = np.array([float(x) for x in fitdata[8].lstrip('[').rstrip(']\n').split(',')])
else:
    Nworld = np.zeros(1, dtype=float)

# Kinematic quantities.
z0 = (1-sqrt(1-t0/tcut))/(1+sqrt(1-t0/tcut))
zQ2 = [(sqrt(tcut+T)-sqrt(tcut-t0))/(sqrt(tcut+T)+sqrt(tcut-t0)) for T in Q2]

# Covariance matrix.
cov = []
for k in range(2*nmax):
    covrow = [float(x) for x in fitdata[9+k].lstrip('[').rstrip(']\n').split(',')]
    cov.append(covrow)
cov = np.array(cov)

# Sum rules.
temp = fitdata[9+2*nmax].replace('[','').replace(']]\n','').split('],')
sumrules = []
for k in range(5):
    sumrules.append([float(x) for x in temp[k].split(',')])
    
# FF values and errors.
GEQ2, GMQ2 = [], []
dGEQ2, dGMQ2 = [], []
for z in zQ2:
    GEQ2.append(array([gecoef_fit[i]*z**i for i in range(len(gecoef_fit))]).sum())
    GMQ2.append(array([gmcoef_fit[i]*z**i for i in range(len(gmcoef_fit))]).sum())
    dGE, dGM = 0, 0
    sr = array([ array([sumrules[0][l+1]*array([k-j for j in range(l)]).prod() for l in range(0,4)]).sum() + array([z**(kmax-n)*(sumrules[4-n][0]*z0**k + array([sumrules[4-n][l+1]*array([k-j for j in range(l)]).prod() for l in range(0,4)]).sum()) for n in reversed(list(range(0,4)))]).sum() for k in range(1,nmax+1)])
    for k in range(nmax):
        for l in range(nmax):
            dGE += cov[k][l]*(z**(k+1)- sumrules[0][0]*z0**(k+1) - sr[k])*(z**(l+1)- sumrules[0][0]*z0**(l+1) - sr[l])
            dGM += cov[nmax+k][nmax+l]*(z**(k+1)- sumrules[0][0]*z0**(k+1) - sr[k])*(z**(l+1)- sumrules[0][0]*z0**(l+1) - sr[l])
    dGE = sqrt(dGE)
    dGM = sqrt(dGM)
    dGEQ2.append(dGE)
    dGMQ2.append(dGM)

GEQ2 = array(GEQ2)
dGEQ2 = array(dGEQ2)
GMQ2 = array(GMQ2)
dGMQ2 = array(dGMQ2)

# Standard dipole form factors and form factor ratios.
Gdip = array([1/(1+T/Lambda2)**2 for T in Q2])
GEdiprat = GEQ2/(GE0*Gdip)
dGEdiprat = dGEQ2/(GE0*Gdip)
GMdiprat = GMQ2/(GM0*Gdip)
dGMdiprat = dGMQ2/(GM0*Gdip)
GEGMrat = (GM0*GEQ2)/(GE0*GMQ2)
dGEGMrat = np.zeros(len(GEGMrat), dtype=float)
for i in range(0, len(GEGMrat)):
    dGEGMrat[i] = GEGMrat[i] * sqrt((dGEQ2[i]/GEQ2[i])**2 + (dGMQ2[i]/GMQ2[i])**2)

# Write output.
of = open(outfile, 'w')
print('##%26s %28s %28s %28s %28s %28s %28s %28s %28s %28s %28s %28s' % ('Q2','z','GEp','dGEp','GMp','dGMp','GEp/GD','dGEp/GD','GMp/muGD','dGMp/muGD','muGE/GM','dGE/GM'), file=of)
for i in range(len(Q2)):
    print('%28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e %28.20e' % (Q2[i], zQ2[i], GEQ2[i], dGEQ2[i], GMQ2[i], dGMQ2[i], GEdiprat[i], dGEdiprat[i], GMdiprat[i], dGMdiprat[i], GEGMrat[i], dGEGMrat[i]), file=of)

print('', file=of)
print('###### Initial Condtion and Fit results #### ', file=of)
print('#1:  Q2, tcut, t0, z0', file=of)
print('#2:  GEp z-expansion parameters', file=of)
print('#3:  GMp z-expansion parameters', file=of)
print('#4:  redchi2, chi2, chi2xsMainz, chi2xsworld, chi2pol, chi2gecoef, chi2gmcoef, chi2csyst, chi2NMainz, chi2Nworld, chi2gefake, chi2gmfake, chi2gefakeHQ, chi2gmfakeHQ, chi2SBS, chi2erad, chi2mrad', file=of)
print('#5:  ndof, num_tot, num_Mainz, num_world, num_pol, num_fake, kmax_temp, kmax_temp, kslope, kMainz, kworld, num_fakeHQ, num_SBS', file=of)
print('#6:  erad, derad, mrad, dmrad', file=of)
print('#7:  Normalization factor for Mainz data', file=of)
print('#8:  Normalization factor for world data', file=of)

print(Q2max, tcut, t0, z0, file=of)
print('  '.join(map(str, gecoef_fit)), file=of)
print('  '.join(map(str, gmcoef_fit)), file=of)
print(redchi2, chi2, chi2xsMainz, chi2xsworld, chi2pol, chi2gecoef, chi2gmcoef, chi2csyst, chi2NMainz, chi2Nworld, chi2gefake, chi2gmfake, chi2gefakeHQ, chi2gmfakeHQ, chi2SBS, chi2erad, chi2mrad, file=of)
print(ndof, num_tot, num_Mainz, num_world, num_pol, num_fake, kmax_temp, kmax_temp, kslope, kMainz, kworld, num_fakeHQ, num_SBS, file=of)
print(erad, derad, mrad, dmrad, file=of)
print('  '.join(map(str, NMainz)), file=of)
print('  '.join(map(str, Nworld)), file=of)
print('', file=of)
print('# 16 x 16 Covariance Matrix', file=of)
print('\n'.join(['  '.join(['{:16}'.format(item) for item in row]) for row in cov]), file=of)

of.close()