############################/
############################/
## Parameterized Form Factor Central Value and Error
############################/
## ID = 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn, 
## Q2 in GeV^2
##
# The parameterization formula returns the uncertainty devided by G(0)*GD, where 
#  GD(Q2) = 1./(1+Q2/0.71)^2
# and GEp(0) =1, GMp(0) = 2.79284356, GEn(0) = 1, GMn(0) = -1.91304272,
#
# The parameterization formula for the Form Factor value is:
#  $$ GN(z) = sum_{i=0}^{N=13}(a_i * z^i) 
#Note that the return value has been divided by (G(Q2=0)*G_Dip)
#
# The parameterization formula for the Form Factor error is:
# $$ log_{10}\frac{\delta G}{G_D} = (L+c_0)\Theta_a(L_1-L) 
#                                 +\sum_{i=1}^{N}(c_i+d_i L)[\Theta_a(L_i-L)-\Theta_a(L_{i+1}-L)]
#                                 +log_{10}(E_{\inf})\Theta_a(L-L_{N+1})$$
# where $L=log_{10}(Q^2)$, $\Theta_{a}(x)=[1+10^{-ax}]^{-1}$. $a=1$.

import numpy as np
from math import *
def GetFF(kID, kQ2):# {{{
    #################################################
    #### z-Expansion Parameters for Form Factor Values
    #################################################{{{
    N=13;#### a0+a[kmax] ,for proton, kmax=12, for neutron, kmax=10 but I set the last two to be zeros
    GEp_Coef_Fit = np.array([0.239235488875,-1.107073054028, 1.452300215920,0.443155941984,-2.372076572589, 1.332410864834, 1.531625962157,-4.256009724873 ,3.818359189884, 1.390170366895,-5.265090474688, 3.621834558198,-0.828842762570])
    GMp_Coef_Fit = np.array([0.264129538390,-1.095990190524, 1.219188836541,0.672800925788,-1.422885032694,-1.412283259774, 1.565906207818, 4.296018791679,-5.616864164282,-2.787434653524, 8.876017255253,-5.889501080558, 1.330896825887])
    GEn_Coef_Fit = np.array([0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0])
    GMn_Coef_Fit = np.array([0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0]) #}}}

    #################################################
    #### Parameters for Form Factor Errors
    #################################################{{{
    Lp= np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
    ## GEp:
    Einf_GEp = 3.157685
    c0_GEp = -0.99959125
    c_GEp = np.array([2.42147739, 0.71594965, 4.15930349, -1.11520973,-0.55258018]) 
    d_GEp = np.array([ -0.57651805, -0.40696275, -1.56967082, 1.06823140, -0.12987175]) 
    ## GMp:
    Einf_GMp = 1.790780
    c0_GMp = -0.69156036
    c_GMp = np.array([0.24782941,2.58825267,1.15520734,5.93012850,3.88750162]) 
    d_GMp = np.array([-1.14687459, -1.33804811,0.65599095, -3.31118121, -1.11984805]) 
    ## GEn:
    Einf_GEn = 6.363856
    c0_GEn =-2.02517380
    c_GEn = np.array([ -0.71438167, 2.17422854, 2.48587860, 2.57114710, 0.56019404])
    d_GEn = np.array([-1.82555135, -0.78957241, -2.86130295, -1.29962824, -0.37809484]) 
    ## GMn:
    Einf_GMn =  4.292528
    c0_GMn = -0.32663855
    c_GMn = np.array([   0.23233237, 5.41338123, 2.51706481, -12.25509149, -5.99225341]) 
    d_GMn = np.array([ -1.16628065, 0.64199162, 6.32905161, 6.22054376 , 1.19648224]) #}}}

    ## Define #{{{
    tcut = 0.0779191396 
    t0 = -0.7 

    Einf = 0.0
    c0 = 0.0
    c = np.zeros(5, dtype=float)
    d = np.zeros(5, dtype=float)
    GN_Coef_Fit = np.zeros(N, dtype=float) #}}}
 
    if kID==1:# {{{
        c0 = c0_GEp
        Einf = Einf_GEp
        c = c_GEp
        d = d_GEp
        GN_Coef_Fit = GEp_Coef_Fit
    elif kID==2:
        c0 = c0_GMp
        Einf = Einf_GMp
        c = c_GMp
        d = d_GMp
        GN_Coef_Fit = GMp_Coef_Fit
    elif kID==3:
        c0 = c0_GEn
        Einf = Einf_GEn
        c = c_GEn
        d = d_GEn
        GN_Coef_Fit = GEn_Coef_Fit
    elif kID==4:
        c0 = c0_GMn
        Einf = Einf_GMn
        c = c_GMn
        d = d_GMn
        GN_Coef_Fit = GMn_Coef_Fit
    else:
        print '*** ERROR, ID is not one of these: 1->GEp, 2->GMp, 3->GEn, 4->GMn'
        GNGD_Fit = -1.0
        GNGD_Err = -2.0
        return GNGD_Fit, GNGD_Err# }}}

    ## Apply the z-expansion formula
    z = (sqrt(tcut+kQ2)-sqrt(tcut-t0))/(sqrt(tcut+kQ2)+sqrt(tcut-t0)) 
    GNQ2 = 0.0
    for i in range(0,N):
        GNQ2 += GN_Coef_Fit[i] * (z**i)
    GDip= 1./(1. + kQ2/0.71)**2
    GNGD_Fit = GNQ2 / GDip #Note that the GN_Coef_Fit has been divided by mu_p or mu_n for GMp and GMn
 
    ## Apply the Parameterization formula for error
    L = log10(kQ2) 
    X0= Lp[0]-L
    THETA0= 1./(1.+10**(-X0))
    TERM1 = (L+c0)*THETA0
    TERM2 = 0.0
    for i in range(0,5):
        X1= Lp[i]-L
        THETA1= 1./(1.+10**(-X1))
        X2= Lp[i+1]-L
        THETA2= 1./(1.+10**(-X2))
        TERM2 += (c[i]+d[i]*L) * ( THETA1 - THETA2 )

    X3= L-Lp[5]
    THETA3= 1./(1.+10**(-X3))
    TERM3 = log10(Einf)*THETA3

    ##LOG10(dG/G(0)/GD)
    GNGD_Err = 10.**(TERM1+TERM2+TERM3)

    return GNGD_Fit, GNGD_Err
# }}}
