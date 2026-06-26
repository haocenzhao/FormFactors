
# coding: utf-8

# # General Parameters

# In[1]:

### Load Python Lib##### 
### plot_FFcompsr_Q2max.py
### Plot form factor ratios vs. Q^2max for fits w/, w/o sum rules with stat errors.
########
import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import FixedLocator, MultipleLocator, FormatStrFormatter
import os, sys, traceback

#very import to enable this so the plots can be showed in the page
get_ipython().magic(u'matplotlib inline')

# Use LaTeX font.
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':20})

import matplotlib.font_manager as font_manager
font_prop = font_manager.FontProperties( size=12)


# In[2]:

## Input Parameters
Q2max = 1000
Q2str = str(Q2max)
Q2plot =10000
kmax = 12
#kmax=int(sys.argv[1])

fitdata = 'world'
#Q2max = int(raw_input('Max Q2 Fit =  '))
#Q2plot = float(raw_input('Max Q2 Plot = '))
#kmax = int(raw_input('Max Z = '))

bnd = '5'
mod = 't0fix7'


# # Load Data

# In[3]:

## Dipole FF
Lambda2 = 0.71
GEn0 = 1
GMn0 = -1.91304272
def GD(Q2):
    return 1./(1+Q2/Lambda2)**2


# ## Load $G_E^n$ Sum-Rules baseline Fitting Results w/ radius constraints

# In[4]:

kmax=10
##Load Q2max=1000 Fitting results#{{{
folder = './z'+str(kmax) # folder containing central fits
#filename=folder+'/Round6_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GEn.dat'
filename=folder+'/out_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GEn.dat'

if os.path.isfile(filename):
    print 'file exist', filename
fitlines = open(filename, 'r').readlines()

##Extract results from text file for non-sum rule fits.
N = 2000
Q2e_n = np.zeros(N, dtype=float)
Ze_n = np.zeros(N, dtype=float)
GE_fit_n = np.zeros(N, dtype=float)
GE_pos_n = np.zeros(N, dtype=float)
GE_neg_n = np.zeros(N, dtype=float)
GErat_fit_n = np.zeros(N, dtype=float)
GErat_pos_n = np.zeros(N, dtype=float)
GErat_neg_n = np.zeros(N, dtype=float)
GED_n = np.zeros(N, dtype=float)
dGEnGD = np.zeros(N, dtype=float)

for i in range(0,N):

    values = fitlines[i+1].split()
    #if float(values[0]) > Q2plot:
    #    break
           
    Q2e_n[i]=float(values[0])
    Ze_n[i]=float(values[1])
    GE_fit_n[i]=float(values[2])
    GE_pos_n[i]=float(values[2])+abs(float(values[3]))
    GE_neg_n[i]=float(values[2])-abs(float(values[3]))
    Gdip = float(values[4])
#    Gdip = 1./(1.+float(values[0])/0.71)**2
    GED_n[i]=Gdip*GEn0
    GErat_fit_n[i]=GE_fit_n[i]/GED_n[i]
    GErat_pos_n[i]=GE_pos_n[i]/GED_n[i]
    GErat_neg_n[i]=GE_neg_n[i]/GED_n[i] 
#}}}


# ## Load $G_E^n$ Sum-Rules  Fitting Results w/o radius constraints

# In[5]:

kmax=10
##Load Q2max=1000 Fitting results#{{{
folder = './z'+str(kmax) # folder containing central fits
#filename=folder+'/Round6_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GEn_noRad.dat'
filename=folder+'/out_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GEn_noRad.dat'

if os.path.isfile(filename):
    print 'file exist', filename
fitlines = open(filename, 'r').readlines()

##Extract results from text file for non-sum rule fits.
N = 2000
Q2e_n1 = np.zeros(N, dtype=float)
Ze_n1 = np.zeros(N, dtype=float)
GE_fit_n1 = np.zeros(N, dtype=float)
GE_pos_n1 = np.zeros(N, dtype=float)
GE_neg_n1 = np.zeros(N, dtype=float)
GErat_fit_n1 = np.zeros(N, dtype=float)
GErat_pos_n1 = np.zeros(N, dtype=float)
GErat_neg_n1 = np.zeros(N, dtype=float)
GED_n1 = np.zeros(N, dtype=float)

for i in range(0,N):

    values = fitlines[i+1].split()
    #if float(values[0]) > Q2plot:
    #    break
           
    Q2e_n1[i]=float(values[0])
    Ze_n1[i]=float(values[1])
    GE_fit_n1[i]=float(values[2])
    GE_pos_n1[i]=float(values[2])+abs(float(values[3]))
    GE_neg_n1[i]=float(values[2])-abs(float(values[3]))
    Gdip = float(values[4])
    #Gdip = 1./(1.+float(values[0])/0.71)**2
    GED_n1[i]=Gdip*GEn0   
    GErat_fit_n1[i]=GE_fit_n1[i]/GED_n1[i]
    GErat_pos_n1[i]=GE_pos_n1[i]/GED_n1[i]
    GErat_neg_n1[i]=GE_neg_n1[i]/GED_n1[i] 

#}}}


# ## Load $G_M^n$ Sum-Rules baseline Fitting Results w/ radius constraints

# In[6]:

kmax=10
##Load Q2max=1000 Fitting results#{{{
folder = './z'+str(kmax) # folder containing central fits
#filename=folder+'/Round6_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GMn.dat'
filename=folder+'/out_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GMn.dat'

if os.path.isfile(filename):
    print 'file exist', filename
fitlines = open(filename, 'r').readlines()

##Extract results from text file for non-sum rule fits.
N = 2000
Q2m_n = np.zeros(N, dtype=float)
Zm_n = np.zeros(N, dtype=float)
GM_fit_n = np.zeros(N, dtype=float)
GM_pos_n = np.zeros(N, dtype=float)
GM_neg_n = np.zeros(N, dtype=float)
GMrat_fit_n = np.zeros(N, dtype=float)
GMrat_pos_n = np.zeros(N, dtype=float)
GMrat_neg_n = np.zeros(N, dtype=float)
GMD_n = np.zeros(N, dtype=float)
dGMnGD = np.zeros(N, dtype=float)

for i in range(0,N):

    values = fitlines[i+1].split()
    #if float(values[0]) > Q2plot:
    #    break
           
    Q2m_n[i]=float(values[0])
    Zm_n[i]=float(values[1])
    GM_fit_n[i]=float(values[2])
    GM_pos_n[i]=float(values[2])+abs(float(values[3]))
    GM_neg_n[i]=float(values[2])-abs(float(values[3]))
    #Gdip = 1./(1.+float(values[0])/0.71)**2
    Gdip = float(values[4])
    GMD_n[i]=Gdip*GMn0       
    GMrat_fit_n[i]=GM_fit_n[i]/GMD_n[i]
    GMrat_pos_n[i]=GM_pos_n[i]/GMD_n[i]
    GMrat_neg_n[i]=GM_neg_n[i]/GMD_n[i] 
#}}}


# ## Load $G_M^n$ Sum-Rules  Fitting Results w/o radius constraints

# In[7]:

kmax=10
##Load Q2max=1000 Fitting results#{{{
folder = './z'+str(kmax) # folder containing central fits
#filename=folder+'/Round6_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GMn_noRad.dat'
filename=folder+'/out_world_sumrules_leastsq_Q2'+Q2str+'_z'+str(kmax)+'_gb'+bnd+'_'+mod+'_GMn_noRad.dat'

if os.path.isfile(filename):
    print 'file exist', filename
fitlines = open(filename, 'r').readlines()

##Extract results from text file for non-sum rule fits.
N = 2000
Q2m_n1 = np.zeros(N, dtype=float)
Zm_n1 = np.zeros(N, dtype=float)
GM_fit_n1 = np.zeros(N, dtype=float)
GM_pos_n1 = np.zeros(N, dtype=float)
GM_neg_n1 = np.zeros(N, dtype=float)
GMrat_fit_n1 = np.zeros(N, dtype=float)
GMrat_pos_n1 = np.zeros(N, dtype=float)
GMrat_neg_n1 = np.zeros(N, dtype=float)
GMD_n1 = np.zeros(N, dtype=float)

for i in range(0,N):

    values = fitlines[i+1].split()
    #if float(values[0]) > Q2plot:
    #    break
           
    Q2m_n1[i]=float(values[0])
    Zm_n1[i]=float(values[1])
    GM_fit_n1[i]=float(values[2])
    GM_pos_n1[i]=float(values[2])+abs(float(values[3]))
    GM_neg_n1[i]=float(values[2])-abs(float(values[3]))
#    Gdip = 1./(1.+float(values[0])/0.71)**2
    Gdip = float(values[4])
    GMD_n1[i]=Gdip*GMn0
    GMrat_fit_n1[i]=GM_fit_n1[i]/GMD_n1[i]
    GMrat_pos_n1[i]=GM_pos_n1[i]/GMD_n1[i]
    GMrat_neg_n1[i]=GM_neg_n1[i]/GMD_n1[i] 

#}}}


# # Calculating Errors 

# In[12]:

## Model Dependent Errors by using different data sets

Q2 = Q2e_n
GE_err_fit = np.zeros(len(Q2e_n), dtype=float)
GM_err_fit = np.zeros(len(Q2e_n), dtype=float)
GE_err_sum = np.zeros(len(Q2e_n), dtype=float)
GM_err_sum = np.zeros(len(Q2e_n), dtype=float)
GE_err_rh = np.zeros(len(Q2e_n), dtype=float)
GM_err_rh = np.zeros(len(Q2e_n), dtype=float)
GE_fit_rh = np.zeros(len(Q2e_n), dtype=float)
GM_fit_rh = np.zeros(len(Q2e_n), dtype=float)

from GetFF import *
for i in range(0, len(Q2e_n)):    
        GE_err_fit[i] = abs(GErat_pos_n[i] - GErat_fit_n[i])# use the central fit values w/ radius constraint
        GE_err_sum[i] = abs(GErat_pos_n[i] - GErat_fit_n[i])# use the error also w/ radius contraints

        GM_err_fit[i] = abs(GMrat_pos_n[i]  - GMrat_fit_n[i]) # use the central fit values w/ radius constraint
        GM_err_sum[i] = abs(GMrat_pos_n1[i] - GMrat_fit_n1[i]) # use the error w/o radius contraints
      
        GE_fit_rh[i], GE_err_rh[i] = GetFF(3, Q2e_n[i])
        GM_fit_rh[i], GM_err_rh[i] = GetFF(4, Q2e_n[i])


# # New Parameterization of Errors

# ## Do some special treatments

# In[13]:

# Fitting the GE Gaps
lnQ2     = np.empty_like(GE_err_sum)
lnGE_Err = np.empty_like(GE_err_sum)
lnGE_DErr = np.empty_like(GE_err_sum)
lnGM_Err = np.empty_like(GM_err_sum)
lnGM_DErr = np.empty_like(GM_err_sum)

lnQ2[0] = log10(Q2[0])
lnGE_Err[0] = log10(GE_err_sum[0])
lnGE_DErr[0] = 0.0212 * lnGE_Err[0] #5%
lnGM_Err[0] = log10(GM_err_sum[0])
lnGM_DErr[0] = 0.0212 * lnGM_Err[0] #5%

for i in range(1, len(GE_err_sum)): 
    lnQ2[i] = log10(Q2[i])
    
    lnGE_Err[i] = log10(GE_err_sum[i])
    lnGE_DErr[i] = 0.0212 #assign 5% in ln(Err) is also 5% coincidently
    lnGE_DErr[i] *= (1.+np.maximum(lnQ2[i]-3.0, 0.0) ) # fixed up to 20GeV^2, 150% or more at 100, 400% or more at 1000
    lnGE_DErr[i] *= lnGE_Err[i]
    
    lnGM_Err[i] = log10(GM_err_sum[i])
    lnGM_DErr[i] = 0.0212 #assign 5% in ln(Err) is also 5% coincidently
    lnGM_DErr[i] *= (1.+np.maximum(lnQ2[i]-3.0, 0.0) ) # fixed up to 20GeV^2, 150% or more at 100, 400% or more at 1000
    lnGM_DErr[i] *= lnGM_Err[i]
    
lnQ2_low = lnQ2[lnQ2<log10(2e-4)]
lnGE_Err_low = lnGE_Err[lnQ2<log10(2e-4)]
lnGE_DErr_low = lnGE_DErr[lnQ2<log10(2e-4)]
lnGM_Err_low = lnGM_Err[lnQ2<log10(2e-4)]
lnGM_DErr_low = lnGM_DErr[lnQ2<log10(2e-4)]

lnQ2_tmp = lnQ2[lnQ2>=log10(2e-4)]
lnGE_Err_tmp = lnGE_Err[lnQ2>=log10(2e-4)]
lnGE_DErr_tmp = lnGE_DErr[lnQ2>=log10(2e-4)]
lnGM_Err_tmp = lnGM_Err[lnQ2>=log10(2e-4)]
lnGM_DErr_tmp = lnGM_DErr[lnQ2>=log10(2e-4)]

lnQ2_cnt = lnQ2_tmp[lnQ2_tmp<=log10(1e3)]
lnGE_Err_cnt = lnGE_Err_tmp[lnQ2_tmp<=log10(1e3)]
lnGE_DErr_cnt = lnGE_DErr_tmp[lnQ2_tmp<=log10(1e3)]
lnGM_Err_cnt = lnGM_Err_tmp[lnQ2_tmp<=log10(1e3)]
lnGM_DErr_cnt = lnGM_DErr_tmp[lnQ2_tmp<=log10(1e3)]


lnQ2_high = lnQ2_tmp[lnQ2_tmp>log10(1e3)]
lnGE_Err_high = lnGE_Err_tmp[lnQ2_tmp>log10(1e3)]
lnGE_DErr_high = lnGE_DErr_tmp[lnQ2_tmp>log10(1e3)]
lnGM_Err_high = lnGM_Err_tmp[lnQ2_tmp>log10(1e3)]
lnGM_DErr_high = lnGM_DErr_tmp[lnQ2_tmp>log10(1e3)]


# ## Check the treatments

# In[15]:

## Check
f1, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,6))
f1.subplots_adjust(bottom=0.08, top=0.96, hspace=0.1)

#######  Subplot for GE/GD.
axes.plot(lnQ2, lnGE_Err, 'k:', color='k',  label='All', linewidth=1.5)
axes.fill_between(lnQ2, lnGE_Err -lnGE_DErr, lnGE_Err+lnGE_DErr, edgecolor='k',facecolor='k',alpha=0.5)

axes.plot(lnQ2_low, lnGE_Err_low, 'k--', color='r',  label='Low-Q2', linewidth=1.5)
axes.fill_between(lnQ2_low, lnGE_Err_low -lnGE_DErr_low, lnGE_Err_low+lnGE_DErr_low, edgecolor='r',facecolor='r',alpha=0.5)

axes.plot(lnQ2_cnt, lnGE_Err_cnt, 'k-', color='b',  label='High-Q2', linewidth=1.5)
axes.fill_between(lnQ2_cnt, lnGE_Err_cnt-lnGE_DErr_cnt, lnGE_Err_cnt+lnGE_DErr_cnt, edgecolor='b',facecolor='b',alpha=0.5)

axes.plot(lnQ2_high, lnGE_Err_high, 'k.', color='g',  label='High-Q2', linewidth=1.5)
axes.fill_between(lnQ2_high, lnGE_Err_high -lnGE_DErr_high, lnGE_Err_high+lnGE_DErr_high, edgecolor='g',facecolor='g',alpha=0.5)

#axes.set_xlim(1e-6, 1e6)
#axes.set_ylim(1e-8, 1e2)


axes.set_xlabel('$ln(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes.set_ylabel('$ln(dG)$', fontsize=30)

#lg= axes.legend(loc='upper left', shadow='true', fontsize='small', numpoints=1)
lg= axes.legend(loc='lower right', shadow='true', fontsize='small', numpoints=1)

axes.xaxis.set_tick_params(width=1, length=7)
axes.yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)


# ## Define Fitting Functions

# In[18]:

##A New Fitting Function
class MyFitter:
    def __init__(self, iter):
        self.iter = iter
    
    def poly1(self, x, par):
        y=np.array([par[j] * x**(j) for j in range(len(par))]).sum()
        return y

    def func(self, params):
        res = []
        chi2=0.0
        for i in range(len(self.X)):
            model = self.poly1(self.X[i],params)
            err = (model-self.Y[i])/self.Err[i]
            chi2 += err*err
            res.append(err)
            self.iter += 1              
        return res
    
    def SetData(self, kX, kY, kErr):
        self.X = kX
        self.Y = kY
        self.Err = kErr
        
    def PrintData(self):
        print 'X = ', self.X
        print 'Y = ',self.Y
    
    def SetIter(self, iter):
        self.iter=iter


# ## Start Fitting low Q2 region ($<2\times 10^{-4}~GeV^2$)

# In[19]:

##### Fitting GE with a polynomial function #############
import scipy.optimize
import time
tol=1e-13
fitmax0 = 2
tot_low = np.ones(fitmax0)

### Now do the fittingfor GE
### Define and initialize parameters
fitter_GE0 = MyFitter(0)
fitter_GE0.SetData(lnQ2_low, lnGE_Err_low, lnGE_DErr_low)

start=time.clock()
GE_bestfit_low, GE_cov_low, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GE0.func, tot_low,full_output=1, ftol=tol)
tot_time = time.clock()-start

GE_chi2_low = (infodict['fvec'][0:len(lnQ2_low)]**2).sum()
GE_ndof_low = len(infodict['fvec'])
print 'For GE low-Q2 Fitting: '
print 'D.O.F = ', GE_ndof_low
print 'Chi2 = ', GE_chi2_low
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg

### Now do the fittingfor GM
### Define and initialize parameters
fitter_GM0 = MyFitter(0)
fitter_GM0.SetData(lnQ2_low, lnGM_Err_low, lnGM_DErr_low)

start=time.clock()
GM_bestfit_low, GM_cov_low, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GM0.func, tot_low,full_output=1, ftol=tol)
tot_time = time.clock()-start

GM_chi2_low = (infodict['fvec'][0:len(lnQ2_low)]**2).sum()
GM_ndof_low = len(infodict['fvec'])
print 'For GM low-Q2 Fitting: '
print 'D.O.F = ', GM_ndof_low
print 'Chi2 = ', GM_chi2_low
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg


# In[20]:

### Read the output parameters
GE_parL = GE_bestfit_low
print GE_parL
lnGE_Err_Fit_low=np.array([fitter_GE0.poly1(lnQ2_low[i], GE_parL) for i in range(len(lnQ2_low))])
GE_Err_Fit_Low = 10**lnGE_Err_Fit_low

GM_parL = GM_bestfit_low
print GM_parL
lnGM_Err_Fit_low=np.array([fitter_GM0.poly1(lnQ2_low[i], GM_parL) for i in range(len(lnQ2_low))])
GM_Err_Fit_Low = 10**lnGM_Err_Fit_low


# ## Start Fitting high Q2 region ($>2\times 10^{-4}~GeV^2$)

# In[21]:

##### Fitting GE with a polynomial function #############
fitmax1 = 2
tot_high = np.ones(fitmax1)

### Now do the fittingfor GE
### Define and initialize parameters
fitter_GE1 = MyFitter(0)
fitter_GE1.SetData(lnQ2_high, lnGE_Err_high, lnGE_DErr_high)

start=time.clock()
GE_bestfit_high, GE_cov_high, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GE1.func, tot_high,full_output=1, ftol=tol)
tot_time = time.clock()-start

GE_chi2_high = (infodict['fvec'][0:len(lnQ2_high)]**2).sum()
GE_ndof_high = len(infodict['fvec'])
print 'For GE high-Q2 Fitting: '
print 'D.O.F = ', GE_ndof_high
print 'Chi2 = ', GE_chi2_high
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg

### Now do the fittingfor GM
### Define and initialize parameters
fitter_GM1 = MyFitter(0)
fitter_GM1.SetData(lnQ2_high, lnGM_Err_high, lnGM_DErr_high)

start=time.clock()
GM_bestfit_high, GM_cov_high, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GM1.func, tot_high,full_output=1, ftol=tol)
tot_time = time.clock()-start

GM_chi2_high = (infodict['fvec'][0:len(lnQ2_high)]**2).sum()
GM_ndof_high = len(infodict['fvec'])
print 'For GM high-Q2 Fitting: '
print 'D.O.F = ', GM_ndof_high
print 'Chi2 = ', GM_chi2_high
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg


# In[22]:

### Read the output parameters
GE_parH = GE_bestfit_high
print GE_parH
lnGE_Err_Fit_high=np.array([fitter_GE1.poly1(lnQ2_high[i], GE_parH) for i in range(len(lnQ2_high))])
GE_Err_Fit_High = 10**lnGE_Err_Fit_high

GM_parH = GM_bestfit_high
print GM_parH
lnGM_Err_Fit_high=np.array([fitter_GM1.poly1(lnQ2_high[i], GM_parH) for i in range(len(lnQ2_high))])
GM_Err_Fit_High = 10**lnGM_Err_Fit_high


# ## Start Fitting central Q2 region ($2\times 10^{-4}\le Q^2 \le 10^{3}~GeV^2$)

# In[23]:

## Replace the first three points with the low-Q2 region points,
## to make sure the fitting go reasonablly at the edge
lnQ2_cnt[0] = lnQ2_low[-1]
lnQ2_cnt[1] = lnQ2_low[-5]
lnQ2_cnt[2] = lnQ2_low[-10]

lnGE_Err_cnt[0] = lnGE_Err_low[-1]
lnGE_Err_cnt[1] = lnGE_Err_low[-5]
lnGE_Err_cnt[2] = lnGE_Err_low[-10]
lnGE_DErr_cnt[0] = 1e-5
lnGE_DErr_cnt[1] = 0.001
lnGE_DErr_cnt[2] = 0.001

lnGM_Err_cnt[0] = lnGM_Err_low[-1]
lnGM_Err_cnt[1] = lnGM_Err_low[-5]
lnGM_Err_cnt[2] = lnGM_Err_low[-10]
lnGM_DErr_cnt[0] = 1e-5
lnGM_DErr_cnt[1] = 0.001
lnGM_DErr_cnt[2] = 0.001

## Replace the first three points with the low-Q2 region points,
## to make sure the fitting go reasonablly at the edge
lnQ2_cnt[-1] = lnQ2_high[1]
lnQ2_cnt[-2] = lnQ2_high[5]
lnQ2_cnt[-3] = lnQ2_high[10]

lnGE_Err_cnt[-1] = lnGE_Err_high[1]
lnGE_Err_cnt[-2] = lnGE_Err_high[5]
lnGE_Err_cnt[-3] = lnGE_Err_high[10]
lnGE_DErr_cnt[-1] = 1e-5
lnGE_DErr_cnt[-2] = 0.001
lnGE_DErr_cnt[-3] = 0.001

lnGM_Err_cnt[-1] = lnGM_Err_high[1]
lnGM_Err_cnt[-2] = lnGM_Err_high[5]
lnGM_Err_cnt[-3] = lnGM_Err_high[10]
lnGM_DErr_cnt[-1] = 1e-5
lnGM_DErr_cnt[-2] = 0.001
lnGM_DErr_cnt[-3] = 0.001


# In[76]:

fitmax2 = 20
tot_cnt = np.ones(fitmax2)


# In[77]:

### First, Add the last 3 points with tiny errors from the low-Q2 fit to make sure the linear behavior at low-Q2
##### Fitting GE with a polynomial function #############
fitter_GE2 = MyFitter(0)
fitter_GE2.SetData(lnQ2_cnt, lnGE_Err_cnt, lnGE_DErr_cnt)

### Now do the fitting
start=time.clock()
GE_bestfit_cnt, GE_cov_cnt, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GE2.func, tot_cnt,full_output=1, ftol=tol)
tot_time = time.clock()-start

GE_chi2_cnt = (infodict['fvec'][0:len(lnQ2_cnt)]**2).sum()
GE_ndof_cnt = len(infodict['fvec'])
print 'For GE high-Q2 Fitting: '
print 'D.O.F = ', GE_ndof_cnt
print 'Chi2 = ', GE_chi2_cnt
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg


##### Fitting GE with a polynomial function #############
fitter_GM2 = MyFitter(0)
fitter_GM2.SetData(lnQ2_cnt, lnGM_Err_cnt, lnGM_DErr_cnt)

### Now do the fitting
start=time.clock()
GM_bestfit_cnt, GM_cov_cnt, infodict, mesg, ier = scipy.optimize.leastsq(fitter_GM2.func, tot_cnt,full_output=1, ftol=tol)
tot_time = time.clock()-start

GM_chi2_cnt = (infodict['fvec'][0:len(lnQ2_cnt)]**2).sum()
GM_ndof_cnt = len(infodict['fvec'])
print 'For GM high-Q2 Fitting: '
print 'D.O.F = ', GM_ndof_cnt
print 'Chi2 = ', GM_chi2_cnt
#print 'Total Time Used:', tot_time
print 'Total Function Called: ', infodict['nfev']
print 'Fit Results: ', ier
print 'Message:', mesg


# In[78]:

### Read the output parameters
GE_parC = GE_bestfit_cnt
print GE_parC
lnGE_Err_Fit_cnt=np.array([fitter_GE2.poly1(lnQ2_cnt[i], GE_parC) for i in range(len(lnQ2_cnt))])

lnGE_Err_Fit_cnt_full=np.array([fitter_GE2.poly1(lnQ2[i], GE_parC) for i in range(len(lnQ2))])
GE_Err_New=np.array([10**(fitter_GE2.poly1(lnQ2[i], GE_parC)) for i in range(len(lnQ2))])

GM_parC = GM_bestfit_cnt
print GM_parC
lnGM_Err_Fit_cnt=np.array([fitter_GM2.poly1(lnQ2_cnt[i], GM_parC) for i in range(len(lnQ2_cnt))])

lnGM_Err_Fit_cnt_full=np.array([fitter_GM2.poly1(lnQ2[i], GM_parC) for i in range(len(lnQ2))])
GM_Err_New=np.array([10**(fitter_GM2.poly1(lnQ2[i], GM_parC)) for i in range(len(lnQ2))])


# ## Save Fitting Results

# In[79]:

## Save parameters

of = open('neutron_fit_par%d.dat'%fitmax2,'w')

print >>of, '======= GE ======'
print >>of, ' Low-Q2: ', GE_parL
print >>of, ' Mid-Q2: ', GE_parC 
print >>of, 'High-Q2: ', GE_parH 
print >>of, ''
print >>of, '======= GM ======'
print >>of, ' Low-Q2: ', GM_parL
print >>of, ' Mid-Q2: ', GM_parC 
print >>of, 'High-Q2: ', GM_parH 
of.close()


# ## Check Fitting Resutls

# In[80]:

f1, axes = plt.subplots(nrows=2, ncols=2, figsize=(14,12))
f1.subplots_adjust(bottom=0.08, top=0.96, hspace=0, wspace=0.3)

#######  Subplot for GE/GD.
axes[0][0].plot(lnQ2, lnGE_Err, 'k:', color='k',  label='All', linewidth=1.5)
#axes[0][0].fill_between(lnQ2, lnGE_Err -lnGE_DErr, lnGE_Err+lnGE_DErr, edgecolor='k',facecolor='k',alpha=0.5)

axes[0][0].plot(lnQ2_low, lnGE_Err_low, 'k--', color='r',  label='Low-Q2', linewidth=1.5)
#axes[0][0].fill_between(lnQ2_low, lnGE_Err_Fit_low -lnGE_DErr_low, lnGE_Err_low+lnGE_DErr_low, edgecolor='r',facecolor='r',alpha=0.5)

axes[0][0].plot(lnQ2_cnt, lnGE_Err_cnt, 'k-', color='b',  label='Mid-Q2', linewidth=1.5)
#axes[0][0].fill_between(lnQ2_cnt, lnGE_Err_Fit_cnt -lnGE_DErr_cnt, lnGE_Err_cnt+lnGE_DErr_cnt, edgecolor='b',facecolor='b',alpha=0.5)

axes[0][0].plot(lnQ2_high, lnGE_Err_high, 'k.', color='g',  label='High-Q2', linewidth=1.5)
#axes[0][0].fill_between(lnQ2_high, lnGE_Err_Fit_high -lnGE_DErr_high, lnGE_Err_high+lnGE_DErr_high, edgecolor='g',facecolor='g',alpha=0.5)

axes[0][0].plot(lnQ2, lnGE_Err_Fit_cnt_full, 'k-', color='skyblue',  label='Mid-Q2 Fit', linewidth=1.5)
axes[0][0].plot(lnQ2_low, lnGE_Err_Fit_low, 'k-', color='maroon',  label='Low-Q2 Fit', linewidth=1.5)
axes[0][0].plot(lnQ2_high, lnGE_Err_Fit_high, 'k-', color='olivedrab',  label='High-Q2 Fit', linewidth=1.5)

axes[0][0].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[0][0].set_ylabel('$log10(dG_E^p/G_D)$', fontsize=30)
axes[0][0].set_xlim(-6, 4)
axes[0][0].set_ylim(-7.9, 3.9)

#lg= axes.legend(loc='upper left', shadow='true', fontsize='small', numpoints=1)
lg= axes[0][0].legend(loc='lower right', shadow='true', fontsize='small', numpoints=1)

axes[0][0].xaxis.set_tick_params(width=1, length=7)
axes[0][0].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

#######  Subplot for Ratio of Fit vs Data
axes[1][0].plot(lnQ2_cnt, lnGE_Err_Fit_cnt/lnGE_Err_cnt, 'k-', color='b',  label='Mid-Q2 Fit (%d Par)'%fitmax2, linewidth=1.5)
axes[1][0].plot(lnQ2_low, lnGE_Err_Fit_low/lnGE_Err_low, 'k-.', color='r',  label='Low-Q2 Fit', linewidth=1.5)
axes[1][0].plot(lnQ2_high, lnGE_Err_Fit_high/lnGE_Err_high, 'k--', color='g',  label='High-Q2 Fit (%d Par)'%fitmax2, linewidth=1.5)

axes[1][0].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[1][0].set_ylabel('$dG^{fit}/dG^{data}$', fontsize=30)
axes[1][0].set_xlim(-6, 4)
axes[1][0].set_ylim(0.8, 1.19)
axes[1][0].grid()

lg1= axes[1][0].legend(loc='lower left', shadow='true', fontsize='small', numpoints=1)
axes[1][0].xaxis.set_tick_params(width=1, length=7)
axes[1][0].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)



#######  Subplot for GE/GD.
axes[0][1].plot(lnQ2, lnGM_Err, 'k-.', color='k',  label='All', linewidth=1.5)
#axes[0][1].fill_between(lnQ2, lnGM_Err -lnGM_DErr, lnGM_Err+lnGM_DErr, edgecolor='k',facecolor='k',alpha=0.5)

axes[0][1].plot(lnQ2_low, lnGM_Err_low, 'k--', color='r',  label='Low-Q2', linewidth=1.5)
#axes[0][1].fill_between(lnQ2_low, lnGM_Err_Fit_low -lnGM_DErr_low, lnGM_Err_low+lnGM_DErr_low, edgecolor='r',facecolor='r',alpha=0.5)


axes[0][1].plot(lnQ2_cnt, lnGM_Err_cnt, 'k-', color='b',  label='Mid-Q2', linewidth=1.5)
#axes[0][1].fill_between(lnQ2_cnt, lnGM_Err_Fit_cnt -lnGM_DErr_cnt, lnGM_Err_cnt+lnGM_DErr_cnt, edgecolor='b',facecolor='b',alpha=0.5)

axes[0][1].plot(lnQ2_high, lnGM_Err_high, 'k.', color='g',  label='High-Q2', linewidth=1.5)
#axes[0][1].fill_between(lnQ2_high, lnGM_Err_Fit_high -lnGM_DErr_high, lnGM_Err_high+lnGM_DErr_high, edgecolor='g',facecolor='g',alpha=0.5)

axes[0][1].plot(lnQ2, lnGM_Err_Fit_cnt_full, 'k-', color='skyblue',  label='Mid-Q2 Fit', linewidth=1.5)
axes[0][1].plot(lnQ2_low, lnGM_Err_Fit_low, 'k-', color='maroon',  label='Low-Q2 Fit', linewidth=1.5)
axes[0][1].plot(lnQ2_high, lnGM_Err_Fit_high, 'k-', color='olivedrab',  label='High-Q2 Fit', linewidth=1.5)

axes[0][1].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[0][1].set_ylabel('$log10(dG_M^p/G_D)$', fontsize=30)
axes[0][1].set_xlim(-6, 4)
axes[0][1].set_ylim(-7.9, 3.9)

#lg= axes.legend(loc='upper left', shadow='true', fontsize='small', numpoints=1)
lg= axes[0][1].legend(loc='lower right', shadow='true', fontsize='small', numpoints=1)

axes[0][1].xaxis.set_tick_params(width=1, length=7)
axes[0][1].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

#######  Subplot for Ratio of Fit vs Data
axes[1][1].plot(lnQ2_cnt, lnGM_Err_Fit_cnt/lnGM_Err_cnt, 'k-', color='b',  label='Mid-Q2 Fit', linewidth=1.5)
axes[1][1].plot(lnQ2_low, lnGM_Err_Fit_low/lnGM_Err_low, 'k-.', color='r',  label='Low-Q2 Fit', linewidth=1.5)
axes[1][1].plot(lnQ2_high, lnGM_Err_Fit_high/lnGM_Err_high, 'k--', color='g',  label='High-Q2 Fit (%d Par)'%fitmax2, linewidth=1.5)

axes[1][1].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[1][1].set_ylabel('$dG^{fit}/dG^{data}$', fontsize=30)
axes[1][1].set_xlim(-6, 4)
axes[1][1].set_ylim(0.8, 1.19)
axes[1][1].grid()

lg1= axes[1][1].legend(loc='lower left', shadow='true', fontsize='small', numpoints=1)
axes[1][1].xaxis.set_tick_params(width=1, length=7)
axes[1][1].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

plt.savefig('GEnGMn_NewPar_%d.pdf'%fitmax2,bbox_inches='tight')


# In[ ]:




# In[ ]:




# In[2]:

import numpy as np
parL = np.zeros((4,2), dtype=float)
parM = np.zeros((4,15), dtype=float)
parH = np.zeros((4,3), dtype=float)
## GEp:
parL[0] = np.array([-1.01098247, 0.99702014]) #Low-Q2
parM[0] = np.array([-2.10178584e+00,   1.33088936e-01,   7.55827806e-01,   7.59738590e-01,
         -3.53382118e-01,  -2.59321009e-01,   2.17394163e-02,  -2.03766951e-02,
        2.09726915e-02,   3.18295827e-02,  -3.21334434e-03,  -7.15373150e-03,
        -4.42359569e-04,   5.39290505e-04,   8.85810466e-05]) #Mid-Q2
parH[0] = np.array([0.74787627, 1.9243734, -0.32761266])   #High-Q2

#GMp:
parL[1] = np.array([-0.71071999, 0.99710867]) #Low-Q2
parM[1] = np.array([-2.04442098e+00,   1.10795466e-01,  -1.27204593e+00,  -4.66266062e-03,
        2.12092278e+00,   2.58935762e-01,  -1.55942082e+00,  -3.17591087e-01,
        5.90266906e-01,   1.81373030e-01,  -1.01765556e-01,  -4.51690195e-02,
        4.00816045e-03,   4.02301079e-03,   4.87577911e-04]) ##Mid-Q2
parH[1] = np.array([0.78050988, 1.97869347, -0.67648286]) #High-Q2

#GEn:
parL[2] = np.array([-2.02311829, 1.00066282]) #Low-Q2
parM[2] = np.array([-2.07343771e+00,   1.13218347e+00,   1.03946682e+00,  -2.79708561e-01,
        -3.39166129e-01,   1.98498974e-01,  -1.45403679e-01,  -1.21705930e-01,
        1.14234312e-01,   5.69989513e-02,  -2.33664051e-02,  -1.35740738e-02,
        7.84044667e-04,   1.19890550e-03,   1.55012141e-04,]) #Mid-Q2:
parH[2] = np.array([0.4553596, 1.95063341, 0.32421279]) #High-Q2:

    #GMn:
parL[3] = np.array([-0.20765505, 0.99767103]) #Low-Q2:
parM[3] = np.array([  -2.07087611e+00,   4.32385770e-02,  -3.28705077e-01,   5.08142662e-01,
        1.89103676e+00,   1.36784324e-01,  -1.47078994e+00,  -3.54336795e-01,
        4.98368396e-01,   1.77178596e-01,  -7.34859451e-02,  -3.72184066e-02,
        1.97024963e-03,   2.88676628e-03,   3.57964735e-04]) #Mid-Q2:
parH[3] = np.array([ 0.50859057, 1.96863291, 0.2321395]) #High-Q2


# In[6]:

for i in range(15):
    print r'& %d &  %12.8e  & %12.8e \\'%(i,parM[2][i], parM[3][i])


# In[ ]:



