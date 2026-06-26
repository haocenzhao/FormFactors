######## 
### FF_sumrules_leastsq.py
### Fits elastic ep-scattering data using scipy.optimize.leastsq.
### Accepts five arguments from command line:
### 1. Q2max for fitting data.
### 2. kmax for order of truncation in z expansion.
### 3. GEbound for bound on parameters of electric form factor.
### 4. Filename for output of full fit results (NB: overwrite).
### 5. Filename to append summarized output.
########
import scipy.optimize # for fitter#{{{
import sys, traceback # for user options from command line
import time
from math import *
from itertools import chain # flatten array of Mainz normalization combinations
import numpy as np
from numpy import append, insert, array, dot
from numpy.linalg import inv # for computing sum rules
from scipy.special import spence # for calculating vertex corrections, NB: spence(x) = Li2(1-x), where Li2 is the dilogarithm/Spence function in (14) of 1307.6227 and in (4.14, B9) of PhysRevC.62.054320
from FF_loaddata import *
from FF_funcs import *#}}}

########################
### CONSTANTS
#########################{{{
pi = 3.14159265358
alpha = 7.29735e-3 # Fine-structure constant.
me = 5.10999e-4 # Mass of the electron in GeV.
Mp = 0.938272 # Mass of the proton in GeV.
Mn = 0.939570 # Mass of the proton in GeV.
mpi = 0.13957 # Mass of the charged pion in GeV.
GE0 = 0. # G_E(0), value of the Sachs electric form factor at q^2 = 0, which for the neutron is 0.
GM0 = -1.91304272 # G_M(0), value of the Sachs magnetic form factor at q^2 = 0, which for the neutron is the magnetic moment.
Lambda2 = 0.71 # Best fit value for parameter in dipole form factor for neutron.
GeVfm = 0.197327 # GeV^(-1) to fm conversion: 1 GeV^(-1) = 0.197327 fm.
barn = 0.389379338e-3 # 1 GeV^(-2)=0.389e-3 barn.

## neutron eletric and magnetic raidd
rE2 = -0.1161 #fm
drE2 = 0.0022 #fm divided by 10 from the original value 0.04 fm
rM = 0.864   #fm^2
drM = 0.0001   #fm^2 divided by 10 from the original value 0.03 fm^2
#########################}}}

########################
### USER OPTIONS
#########################{{{
# Q2max
Q2max = float(sys.argv[1])
# kmax of fit, nmax is free parameters after 4 sum rules are imposed
kmax = int(sys.argv[2])
nmax = kmax-4
# option for datasets to fit: 1 GEn, 2 GMn
gnopt = int(sys.argv[3])
# option for TPE correction: 1 for Feshbach, 2 for none, default is Blunden.
tpeopt = 0#int(sys.argv[4])
# Kinematic variable for systematics using lambda function, with input one Mainz data point d, an array of the form [beam energy, Q2, tau, theta, z, cs, dcs, norms, spec].

t0mode=sys.argv[4]

# tcut for z-expansion
tcut = 4*mpi**2
# t0 and z0 for z-expansion
t0 = 0
if t0mode=='_t0zero': 
    t0 = 0
elif t0mode=='_t0opt': 
    t0 = tcut*(1-sqrt(1+Q2max/tcut))
elif t0mode=='_t0fix': 
    t0 = -0.8
elif t0mode=='_t0fix7': 
    t0 = -0.7
elif t0mode=='_t0fix4': 
    t0 = -0.4
else:
    print('t0 modle is unknown (%s), use t0=0'%t0mode)

z0 = (1-sqrt(1-t0/tcut))/(1+sqrt(1-t0/tcut))

# bounds
GEbound = int(sys.argv[5])
GMbound = GEbound*GM0

kinvar = lambda d : d[3]
# tol parameter for chi^2 (see scipy.optimize.leastsq definition of ftol)
tol = 1e-14

GNbound =0.0
if gnopt==1:
    GN0=GE0
    GNbound = GEbound
elif gnopt==2:
    GN0=GM0
    GNbound = GMbound
else:
    print("*** ERROR, wrong data option (1 or 2) =", gnopt)

fakeopt = 1 ##add fake points to contrain radii

rE_opt = 1
rM_opt = 1
if( 'noRad' in sys.argv[6]):
    rE_opt = 0
    rM_opt = 0

if rE_opt==1:
    print('---- using rE radius constraint, rE2_fixed = %4.3f+/-%4.3f'%(rE2, drE2)) 
else:
    print('---- No rE constraint! ')
if rM_opt==1:
    print('---- using rM radius constraint, rM_fixed = %4.3f+/-%4.3f'%(rM, drM)) 
else:
    print('---- No rM constraint! ')
#########################}}}

########################
### Use matrix of sum rules, srmat, to compute a0, a_{kmax-n} for n = {0,1,2,3}
### for a form factor using the z expansion, with given z0 and nmax = kmax - 4.
### GE0 = GE(Q2=0), parg = [a_1, ..., a_nmax].
########################{{{
def applysumrules(srmat, z0, nmax, GN0, parg):

    # Compute combinations of parameters for sum rules.
    vG = [GN0, 0, 0, 0, 0]
    zn0 = z0
    for i in range(nmax):
        vG[0] -= parg[i]*zn0
        vG[1] -= parg[i]
        vG[2] -= (i+1)*parg[i]
        vG[3] -= (i+1)*i*parg[i]
        vG[4] -= (i+1)*i*(i-1)*parg[i]
        zn0 = z0*zn0

    # z expansion sum rules: compute a0, a_{kmax-n} for n = {0,1,2,3}.
    a0 = sum(srmat[0][j]*vG[j] for j in range(len(vG)))
    for i in range(1,len(srmat)):
        parg = append(parg, sum(srmat[i][j]*vG[j] for j in range(len(vG))))
    
    return a0, parg
########################}}}


########################
### Fitting function computing residuals for data points and Gaussian bounds.
### Assume params = [a_1, ..., a_nmax, N_1, ..., N_GN, where
### a_i are GE coefficients,
### and N's are normalization parameters for world data.
#########################{{{
def fitfuncsumrules(params):
    # Iteration number and array for residuals.
    global iter
    res = []

    count1= 0
    count2= 0
    count3= 0
    count4= 0
    count5= 0

    # Extract parameters.
    gncoef = params[0:nmax]
    Ngn = params[nmax:]
    print(nmax, params)

    # Apply sum rules.
    a0, gncoef = applysumrules(srmat, z0, nmax, GN0, gncoef)
    
    # GEn or GMn residuals.
    chi2gn = 0
    for d in dataGN:
        # Multiply FF and dFF by normalization parameter.
        norm = 1.0
        if len(normGN)>0 and d[5]>1e-5:
            norm = Ngn[normGN.index([d[4],d[5]])]
            #print 'NORM=',Ngn, d[4], normGN.index(d[4]), norm
        data_gn = d[2]*norm
        datae_gn = d[3]*norm

        # Calculate model FF value.
        model = getff(d[0], d[1], a0, gncoef)
        if gnopt==2:
            gd = (1/(1+d[0]/Lambda2)**2) # dipole form factor.
            model = model/gd/GM0 # GMn data was normalized to GD

        # Compute residuals and FF chi^2.
        err = (data_gn - model)/datae_gn
        chi2gn += err*err
        res.append(err)
        count1+=1

    # Normalization residuals.
    chi2Ngn = 0
    for k in range(len(Ngn)):
    #for k in range(len(normGN)):
        #err = (Ngn[k]-1)/0.02  ## make a 2% normalzation residual, 08/08/2016
        err = 0.0
        if normGN[k][1]>1e-33:#remove 0 normalization value
            err = (Ngn[k]-1)/normGN[k][1]
        res.append(err)
        chi2Ngn += err*err
        count2+=1

##    Form factor coefficient residuals for [a1, ... anmax]
    chi2gncoef = 0
    if abs(GNbound) > 0:
        for a in gncoef[0:kmax]:#[0:nmax]:
            err = a/GNbound
            res.append(err)
            count3+=1
            chi2gncoef += err*err

# fake GEn or GMn point for high-Q2 contraint
    chi2gnfake = 0
    for d in datafake:
        err = 0.0
        model = getff(d[0], d[1], a0, gncoef);
        if gnopt==2:
            gd = (1/(1+d[0]/Lambda2)**2) # dipole form factor.
            model = model/gd/GM0 # GMn data was normalized to GD
      # Compute residuals and fake chi^2.
        err= (model - d[2])/d[3]
        chi2gnfake += err*err
        res.append(err)
        count4+=1
        
# SBS GEn points
    chi2SBS_GEn = 0
    for d in dataSBS_GEn:
        model = getff(d[0], d[1], a0, gncoef)
        err = (model - d[2]) / d[3]
        chi2SBS_GEn += err * err
        res.append(err)
        count5 += 1

# radius contraint with the exisiting results
    chi2rad = 0
    radius2= zrad2(t0, tcut, gncoef) #return the R^2 value
    err = 0.0
    if gnopt==1 and rE_opt==1:
        err = (radius2- rE2)/drE2 ## actually r^2 for GEn
    if gnopt==2 and rM_opt==1:
        radius2/=GN0
        err = (sqrt(radius2) - rM)/drM
    chi2rad += err*err
    res.append(err)

    # Print fit update every 50 calls.
    if iter %5 == 0:
        print(iter, 'chi2gn=%f, chi2Ngn=%f, chi2R=%f, chi2gncoef=%f, chi2SBS=%f, dof=%d, %d, %d, %d, %d' %(chi2gn, chi2Ngn, chi2rad, chi2gncoef, chi2SBS_GEn, count1, count2, count3, count4, count5))

    iter += 1
    return res
#########################}}}

########################
### EXTRACT DATA
#########################{{{
# Extract GEn or GMn fit data.
dataGN, expGN, normGN = loadGN(Q2max, t0, tcut, gnopt)
numGN, kGN, kNorm = len(dataGN), len(expGN), len(normGN)
print("Loaded", numGN, "GEn or GMn data points from", kGN,'experiments, with',kNorm,' normalization!')
# print "normGN=", normGN

#temperately disable the normalization fitting for GEn
if gnopt==1:
    normGN=[]
    kNorm=0

# Extract fake data.
if fakeopt == 1:
    #datafake= loadgnfake(Q2max, t0, tcut, gnopt)
    datafake= loadgnfakeHQ(Q2max, t0, tcut, gnopt)
    numfake = len(datafake)
    print(datafake)
else:
    datafake = []
    numfake = 0
if( 'noHQ' in sys.argv[6]):
    datafake = []
    numfake = 0
print("Loaded", numfake, "fake data points.")

# Extract SBS GEn data.
if gnopt == 1:
    dataSBS_GEn = load_SBS_GEn(Q2max, t0, tcut)
    numSBS_GEn = len(dataSBS_GEn)
else:
    dataSBS_GEn = []
    numSBS_GEn = 0
print("Loaded", numSBS_GEn, "SBS GEn data points.")

#}}}

########################
### PERFORM FIT
#########################{{{
# Compute sum rules matrix.
#srmat = sumrules(z0, nmax) # sum rules matrix
srmat = inv([[1,z0**(nmax+1),z0**(nmax+2),z0**(nmax+3),z0**(nmax+4)], [1,1,1,1,1], [0,nmax+1,nmax+2,nmax+3,nmax+4], [0,(nmax+1)*nmax,(nmax+2)*(nmax+1),(nmax+3)*(nmax+2),(nmax+4)*(nmax+3)], [0,(nmax+1)*nmax*(nmax-1),(nmax+2)*(nmax+1)*nmax,(nmax+3)*(nmax+2)*(nmax+1),(nmax+4)*(nmax+3)*(nmax+2)]])

method = int(sys.argv[7])
print("method = ", method)
if kmax<=6:
    ###initial values for GEn fit
    if t0mode=='_t0zero' and gnopt==1: 
        # z-exp seed, t0 = 0
        gncoef_init = [ 0.00,0.08] + [0]*(nmax-2)
    if ('_t0opt' in t0mode or '_t0fix' in t0mode) and gnopt==1: 
        # z-exp seed, t0 = t0opt
        gncoef_init = [ 0.045,0.035] + [0]*(nmax-2)

    ###initial values for GMn fit
    if t0mode=='_t0zero' and gnopt==2: 
        # z-exp seed, t0 = 0
        gncoef_init = [-1.91,1.66] + [0]*(nmax-2)
    if ('_t0opt' in t0mode or '_t0fix' in t0mode) and gnopt==2: 
        # z-exp seed, t0 = t0opt
        gncoef_init = [-0.43,1.97] + [0]*(nmax-2)
else:
    prv_file= sys.argv[6]
    prv_file= prv_file.replace('Apr2026_world_','jun14_world_')
    prv_file= prv_file.replace('Apr2026_world_','jun14_world_')
    prv_file= prv_file.replace('Q2%d_'%Q2max,'Q21000_')
    if( 'noHQ' in sys.argv[6]):
        prv_file= prv_file.replace('_noHQ.dat','.dat')
    if( 'noRad' in sys.argv[6]):
        prv_file= prv_file.replace('_noRad.dat','.dat')
    if( '_t0opt' in sys.argv[6]):
        prv_file= prv_file.replace('_t0opt','_t0fix7')
    if( '_t0zero' in sys.argv[6]):
        prv_file= prv_file.replace('_t0zero','_t0fix7')

    if method ==1:
        ##method-1: read the fitting parameters from previous fit
        prv_file = prv_file.replace('_z%d_gb'%(kmax),'_z%d_gb'%(kmax-1))
        prv_file = prv_file.replace('/z%d/'%(kmax),'/z%d/'%(kmax-1))
    elif method ==2:
        ##method-2: read the fitting parameters from t0opt fit for t0zero fit
        prv_file = prv_file.replace('_opt','_t0fix7')
    elif method ==3:
        ##method-3: read the fitting parameters from bound-fit for unbound-fit
        prv_file = prv_file.replace('_bound0','_bound5')
        prv_file = prv_file.replace('_gb0','_gb5')
    else:
        ##method-1 as the default
        prv_file = prv_file.replace('_z%d_gb'%(kmax),'_z%d_gb'%(kmax-1))
        prv_file = prv_file.replace('/z%d/'%(kmax),'/z%d/'%(kmax-1))

    print('-- Getting parameters from old_file = ', prv_file)
    prv_fit= open(prv_file, 'r').readlines()
    gncoef_temp = [float(x) for x in prv_fit[3].lstrip('[').rstrip(']\n').split(',')]
    if method ==1:
        gncoef_init = [gncoef_temp[i] for i in range(1, len(gncoef_temp)-4) ] + [0]*1
    else:
        gncoef_init = [gncoef_temp[i] for i in range(1, len(gncoef_temp)-4) ]


iter = 0 # number of fit iterations
start = time.perf_counter() # time for diagnostics
result, cov, infodict, mesg, ier = scipy.optimize.leastsq(fitfuncsumrules, gncoef_init+[1]*kNorm, full_output=1, ftol=tol) # can set ftol, xtol parameters for tolerance in chi^2 and parameter values
elapsed = time.perf_counter() - start
#}}}

########################
### EXTRACTION AFTER FITTING
#########################{{{
# Parameter arrays for fit results.
gncoef = [] 
Ngn = []

# Extract form factor coefficients; compute a_0 and rE
for i in range(nmax):
    gncoef.append(result[i])
a0, gncoef = applysumrules(srmat, z0, nmax, GN0, gncoef)

gncoef = list(gncoef)

# Extract correlated systematic slopes and normalizations.
for i in range(kNorm):
    Ngn.append(result[nmax+i])

# Extract chi-square values.
istart, iend = 0, numGN
chi2Num = (infodict['fvec'][istart:iend]**2).sum()
istart, iend = iend, iend+kNorm
chi2Ngn = (infodict['fvec'][istart:iend]**2).sum()
istart, iend = iend, iend+kmax
chi2gncoef = (infodict['fvec'][istart:iend]**2).sum()
istart, iend = iend, iend+numfake
chi2gnfake = (infodict['fvec'][istart:iend]**2).sum()
istart, iend = iend, iend+numSBS_GEn
chi2SBS_GEn = (infodict['fvec'][istart:iend]**2).sum()
istart, iend = iend, iend+1
chi2rad = (infodict['fvec'][istart:iend]**2).sum()

chi2 = chi2Num + chi2Ngn + chi2gncoef + chi2gnfake + chi2SBS_GEn + chi2rad
print('---Len', iend, len(infodict['fvec'][:]**2))
#}}}

########################
### RADIUS ERRORS FROM COVARIANCE MATRIX
#########################{{{
# Follow Section 15.6, especially Eq. 15.6.5, of Numerical Recipes applied to <r^2>.
sr = array([ array([(kmax-n)*z0**(kmax-1-n)*(srmat[4-n][0]*z0**k + array([srmat[4-n][l+1]*array([k-j for j in range(l)]).prod() for l in range(0,4)]).sum() ) for n in reversed(list(range(0,4)))]).sum() for k in range(1,nmax+1)])
cov = array([[cov[i][j] for j in range(nmax)] for i in range(nmax)])

rad2= 0.0
drad2= 0.0
if gnopt==1:
    r2lindep = array([k*z0**(k-1) - sr[k-1] for k in range(1,nmax+1)])*(-6)/tcut*sqrt(1-t0/tcut)/(1+sqrt(1-t0/tcut))**2*GeVfm**2
    r2lindep = r2lindep
    rad2= zrad2(t0, tcut, gncoef) #return the R^2 value
    if abs(rad2)>1e-33:
        drad2= sqrt(dot(r2lindep, dot(cov, r2lindep)))
    else:
        drad2 = 0.0
elif gnopt==2:
    r2lindep = array([k*z0**(k-1) - sr[k-1] for k in range(1,nmax+1)])*(-6)/tcut*sqrt(1-t0/tcut)/(1+sqrt(1-t0/tcut))**2*GeVfm**2
    r2lindep = r2lindep/GN0
    rad2= zrad2(t0, tcut, gncoef)/GN0  #return the R^2 value
    if abs(rad2)>1e-33:
        drad2=     (dot(r2lindep, dot(cov, r2lindep)))/4/rad2
        # drad2= sqrt(dot(r2lindep, dot(cov, r2lindep)))
    else:
        drad2 = 0.0
else:
    rad2= 0.0
    drad2= 0.0

if np.isinf(rad2):
    rad2= 0.0
if np.isinf(drad2):
    drad2= 0.0
#}}}

##make sure to insert a0 after the radius calculation is done
gncoef.insert(0, a0) # [a_0, a_1, ..., a_kmax]

########################
### OUTPUT TO SCREEN
#########################{{{
# Print fit attributes.
print('Time for fit:', elapsed)
print('Number of function calls:', infodict['nfev'])
print('scipy.optimize fit result:', ier)
print('scipy.optimize fit message:', mesg)

# Print chi^2 information.
print('Chi^2 tot:', chi2)
print('Chi^2 Norm:', chi2Ngn)
print('Chi^2 Bound:', chi2gncoef)
print('Chi^2 fake', chi2gnfake)
print('Chi^2 SBS GEn', chi2SBS_GEn)
print('Chi^2 Rad', chi2rad)
print('Total no. data points:', numGN + numSBS_GEn)
print('Total no. fit parameters:', nmax)
print('Total no. experiments:', kGN)
print('Total no. normalizations:', kNorm)
print('No. SBS GEn points:', numSBS_GEn)
ndof = len(infodict['fvec'])-nmax
print('Total degrees of freedom:', ndof)
redchi2 = chi2/ndof
print('Reduced chi^2:', redchi2)

# Print coefficients and normalizations.
print('a_i:', gncoef)

if gnopt==1:
    print('r2_E, dr2_E :', rad2, drad2)
if gnopt==2:
    print('rM, d_M :', sqrt(rad2), sqrt(drad2))
    #print 'rM, d_M :', sqrt(rad2), drad2/2/sqrt(rad2)
#}}}

########################
### OUTPUT TO FILE
#########################{{{
of = open(sys.argv[6],"w")
print(Q2max, tcut, t0, GEbound, GMbound, gnopt, tpeopt, tol, file=of)
print(ndof, nmax, numGN, kGN, kNorm, numSBS_GEn, file=of)
print(redchi2, chi2, chi2Ngn, chi2gncoef, chi2gnfake, chi2SBS_GEn, chi2rad, file=of)
print(gncoef, file=of)
print(rad2, drad2, file=of)

# Covariance matrix and sum rules
for i in range(nmax):
    covrow = []
    for j in range(nmax):
        covrow.append(cov[i][j]) #!May 23, *redchi2) # need to multiply cov by reduced chi^2 to get true approximation to Hessian
    print(covrow, file=of)
print(list(list(srmat[i]) for i in range(len(srmat))), file=of)

# Fit attributes, initializations, normalizations
print(infodict['nfev'], elapsed, ier, file=of)
print(gncoef_init, file=of)
#}}}
of.close()