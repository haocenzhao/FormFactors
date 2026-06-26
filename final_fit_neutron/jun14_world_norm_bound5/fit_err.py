### Load Python Lib##### # {{{
########
import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import FixedLocator, MultipleLocator, FormatStrFormatter
import os, sys, traceback

#very import to enable this so the plots can be showed in the page
# get_ipython().magic(u'matplotlib inline')

# Use LaTeX font.
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':20})

import matplotlib.font_manager as font_manager
font_prop = font_manager.FontProperties( size=12)# }}}

target = sys.argv[1]
# # Load Lookup Tables for Errors # {{{
filename='%s_baseline_seg5.dat'%target

if os.path.isfile(filename):
    print 'file exist', filename

fitlines = open(filename, 'r').readlines()

##Extract results from text file for non-sum rule fits.
N = 2000
Q2 = np.zeros(N, dtype=float)
GE_fit_rat = np.zeros(N, dtype=float)
GM_fit_rat = np.zeros(N, dtype=float)
GE_err_sum = np.zeros(N, dtype=float)
GM_err_sum = np.zeros(N, dtype=float)
dGEGDn = np.zeros(N, dtype=float)
dGMGDn = np.zeros(N, dtype=float)

for i in range(0,N):
    values = fitlines[i+1].split()
           
    Q2[i]=float(values[0])
    GE_fit_rat[i]=float(values[1])
    GE_err_sum[i] = abs(float(values[2]) ) #/ abs(float(values[2]))
    dGEGDn[i] = abs(float(values[3]) )
    GM_fit_rat[i]=float(values[4])
    GM_err_sum[i] = abs(float(values[5]) ) #/ abs(float(values[4]))
    dGMGDn[i] = abs(float(values[6]) )  
#}}}
       
######################################
# # New Parameterization of Errors
######################################

# ## Do some special treatments# {{{
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

GE_err_new = GE_err_sum.copy()
for i in range(1, len(GE_err_sum)): 
    lnQ2[i] = log10(Q2[i])

    if target=='proton':#Smooth out the dips of GEp
        GE_err_new[i] = np.maximum(GE_err_new[i-1],GE_err_new[i]) 
    
    lnGE_Err[i] = log10(GE_err_new[i])
    lnGE_DErr[i] = 0.0212 #assign 5% in ln(Err) is also 5% coincidently
    lnGE_DErr[i] *= (1.+np.maximum(lnQ2[i]-3.0, 0.0) ) # fixed up to 20GeV^2, 150% or more at 100, 400% or more at 1000
    lnGE_DErr[i] *= lnGE_Err[i]
    
    lnGM_Err[i] = log10(GM_err_sum[i])
    lnGM_DErr[i] = 0.0212 #assign 5% in ln(Err) is also 5% coincidently
    lnGM_DErr[i] *= (1.+np.maximum(lnQ2[i]-3.0, 0.0) ) # fixed up to 20GeV^2, 150% or more at 100, 400% or more at 1000
    lnGM_DErr[i] *= lnGM_Err[i]
    
lnQ2_low = lnQ2[lnQ2<log10(1e-4)]
lnGE_Err_low = lnGE_Err[lnQ2<log10(1e-4)]
lnGE_DErr_low = lnGE_DErr[lnQ2<log10(1e-4)]
lnGM_Err_low = lnGM_Err[lnQ2<log10(1e-4)]
lnGM_DErr_low = lnGM_DErr[lnQ2<log10(1e-4)]

lnQ2_tmp = lnQ2[lnQ2>=log10(1e-4)]
lnGE_Err_tmp = lnGE_Err[lnQ2>=log10(1e-4)]
lnGE_DErr_tmp = lnGE_DErr[lnQ2>=log10(1e-4)]
lnGM_Err_tmp = lnGM_Err[lnQ2>=log10(1e-4)]
lnGM_DErr_tmp = lnGM_DErr[lnQ2>=log10(1e-4)]

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
lnGE_DErr_high[0] = 1e-10
lnGM_DErr_high[0] = 1e-10# }}}

# ## Define Fitting Functions# {{{
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
# }}}

######################################
# ## Start Fitting low Q2 region ($<1\times 10^{-4}~GeV^2$)
######################################

##### Fitting GE&GM with a polynomial function ############## {{{
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

# }}}

### Read the output parameters# {{{
GE_parL = GE_bestfit_low
print GE_parL
lnGE_Err_Fit_low=np.array([fitter_GE0.poly1(lnQ2_low[i], GE_parL) for i in range(len(lnQ2_low))])
GE_Err_Fit_Low = 10**lnGE_Err_Fit_low

GM_parL = GM_bestfit_low
print GM_parL
lnGM_Err_Fit_low=np.array([fitter_GM0.poly1(lnQ2_low[i], GM_parL) for i in range(len(lnQ2_low))])
GM_Err_Fit_Low = 10**lnGM_Err_Fit_low
# }}}

######################################
# ## Start Fitting high Q2 region ($>2\times 10^{-4}~GeV^2$)
######################################

##### Fitting GE&GM with a polynomial function ############## {{{
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
# }}}

### Read the output parameters# {{{
GE_parH = GE_bestfit_high
print GE_parH
lnGE_Err_Fit_high=np.array([fitter_GE1.poly1(lnQ2_high[i], GE_parH) for i in range(len(lnQ2_high))])
GE_Err_Fit_High = 10**lnGE_Err_Fit_high

GM_parH = GM_bestfit_high
print GM_parH
lnGM_Err_Fit_high=np.array([fitter_GM1.poly1(lnQ2_high[i], GM_parH) for i in range(len(lnQ2_high))])
GM_Err_Fit_High = 10**lnGM_Err_Fit_high
# }}}

######################################
# ## Start Fitting central Q2 region ($2\times 10^{-4}\le Q^2 \le 10^{3}~GeV^2$)
######################################

## Replace the first three points with the low-Q2 region points,# {{{
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
lnGE_DErr_cnt[-1] = 1e-5*lnGE_Err_high[1]
lnGE_DErr_cnt[-2] = 1e-5*lnGE_Err_high[5]
lnGE_DErr_cnt[-3] = 1e-5*lnGE_Err_high[10]

lnGM_Err_cnt[-1] = lnGM_Err_high[1]
lnGM_Err_cnt[-2] = lnGM_Err_high[5]
lnGM_Err_cnt[-3] = lnGM_Err_high[10]
lnGM_DErr_cnt[-1] = 1e-5*lnGM_Err_high[1]
lnGM_DErr_cnt[-2] = 1e-5*lnGM_Err_high[5]
lnGM_DErr_cnt[-3] = 1e-5*lnGM_Err_high[10]# }}}

fitmax2 = int(sys.argv[1])
tot_cnt = np.ones(fitmax2)

##### Fitting GE with a polynomial function ############## {{{
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
print 'Message:', mesg# }}}

##### Fitting GM with a polynomial function ############## {{{
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
print 'Message:', mesg# }}}

### Read the output parameters# {{{
GE_parC = GE_bestfit_cnt
print GE_parC
lnGE_Err_Fit_cnt=np.array([fitter_GE2.poly1(lnQ2_cnt[i], GE_parC) for i in range(len(lnQ2_cnt))])

lnGE_Err_Fit_cnt_full=np.array([fitter_GE2.poly1(lnQ2[i], GE_parC) for i in range(len(lnQ2))])
GE_Err_New=np.array([10**(fitter_GE2.poly1(lnQ2[i], GE_parC)) for i in range(len(lnQ2))])

GM_parC = GM_bestfit_cnt
print GM_parC
lnGM_Err_Fit_cnt=np.array([fitter_GM2.poly1(lnQ2_cnt[i], GM_parC) for i in range(len(lnQ2_cnt))])

lnGM_Err_Fit_cnt_full=np.array([fitter_GM2.poly1(lnQ2[i], GM_parC) for i in range(len(lnQ2))])
GM_Err_New=np.array([10**(fitter_GM2.poly1(lnQ2[i], GM_parC)) for i in range(len(lnQ2))])# }}}

######################################
# ## Save Fitting Results
######################################
## Save parameters# {{{
of = open('%s_fit_par%d.dat'%(target,fitmax2),'w')

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
# }}}

######################################
# ## Check Fitting Resutls
######################################
f1, axes = plt.subplots(nrows=2, ncols=2, figsize=(14,12))# {{{
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
axes[0][0].set_ylabel('$log10(dG_E^%s/G_D)$'%TGT, fontsize=30)
axes[0][0].set_xlim(-6, 4)
axes[0][0].set_ylim(-7.9, 3.9)

#lg= axes.legend(loc='upper left', shadow='true', fontsize='small', numpoints=1)
lg= axes[0][0].legend(loc='lower right', shadow='true', fontsize='small', numpoints=1)

axes[0][0].xaxis.set_tick_params(width=1, length=7)
axes[0][0].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

#######  Subplot for Ratio of Fit vs Data
axes[1][0].plot(lnQ2_cnt, 10**lnGE_Err_Fit_cnt/10**lnGE_Err_cnt, 'k-', color='b',  label='Mid-Q2 Fit (%d Par)'%fitmax2, linewidth=1.5)
axes[1][0].plot(lnQ2_low, 10**lnGE_Err_Fit_low/10**lnGE_Err_low, 'k-.', color='r',  label='Low-Q2 Fit (2 Par)', linewidth=1.5)
axes[1][0].plot(lnQ2_high, 10**lnGE_Err_Fit_high/10**lnGE_Err_high, 'k--', color='g',  label='High-Q2 Fit (2 Par)', linewidth=1.5)

axes[1][0].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[1][0].set_ylabel('$dG^%s_E(fit)/dG^%s_E(data)$'%(TGT,TGT), fontsize=30)
axes[1][0].set_xlim(-6, 4)
axes[1][0].set_ylim(0.8, 1.20)
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
axes[0][1].set_ylabel('$log10(dG_M^%s/G_D)$'%TGT, fontsize=30)
axes[0][1].set_xlim(-6, 4)
axes[0][1].set_ylim(-7.9, 3.9)

#lg= axes.legend(loc='upper left', shadow='true', fontsize='small', numpoints=1)
lg= axes[0][1].legend(loc='lower right', shadow='true', fontsize='small', numpoints=1)

axes[0][1].xaxis.set_tick_params(width=1, length=7)
axes[0][1].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

#######  Subplot for Ratio of Fit vs Data
axes[1][1].plot(lnQ2_cnt, 10**lnGM_Err_Fit_cnt/10**lnGM_Err_cnt, 'k-', color='b',  label='Mid-Q2 Fit (%d Par)'%fitmax2, linewidth=1.5)
axes[1][1].plot(lnQ2_low, 10**lnGM_Err_Fit_low/10**lnGM_Err_low, 'k-.', color='r',  label='Low-Q2 Fit (2 Par)', linewidth=1.5)
axes[1][1].plot(lnQ2_high, 10**lnGM_Err_Fit_high/10**lnGM_Err_high, 'k--', color='g',  label='High-Q2 Fit (2 Par)', linewidth=1.5)

axes[1][1].set_xlabel('$log10(Q^2)$ $[\mathrm{GeV}^2]$', fontsize=30)
axes[1][1].set_ylabel('$dG^%s_M(fit)/dG^%s_M(data)$'%(TGT,TGT), fontsize=30)
axes[1][1].set_xlim(-6, 4)
axes[1][1].set_ylim(0.8, 1.20)
axes[1][1].grid()

lg1= axes[1][1].legend(loc='lower left', shadow='true', fontsize='small', numpoints=1)
axes[1][1].xaxis.set_tick_params(width=1, length=7)
axes[1][1].yaxis.set_tick_params(width=1, length=7)
plt.tick_params(labelsize=30)

plt.savefig('%s_NewPar_%d.pdf'%(target,fitmax2),bbox_inches='tight')
# }}}
