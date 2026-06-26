######## 
### FF_funcs.py
### Contains functions to compute:
### 1. Various kinematic quantities in electron-proton (or generic heavy particle) scattering.
### 2. Form factors for z-expansion and (inverse) polynomial functional forms and their radii.
### 3. Form factor ratios and reduced cross sections.
### 4. Sum rules for z expansion.
########

from math import *
from numpy import append, insert, array, dot

########################
### CONSTANTS
########################
GeVfm = 0.197327 # GeV^(-1) to fm conversion: 1 GeV^(-1) = 0.197327 fm.
########################

########################
### 1. KINEMATIC QUANTITIES
########################

########################
### Q2 for fixed value of energy E, angle th, heavy particle mass M.
########################
def Qsq(E, th, M):
    return 2*E**2/(1/(1-cos(th)) + E/M)
########################

########################
### Tau for fixed value of Q2, heavy particle mass M.
########################
def tau(Q2, M):
    return Q2/4./M**2;
########################

########################
### Final electron scattering angle for fixed value of Q2, initial energy E, heavy particle mass M.
########################
def theta(E, Q2, M):
    return 2*acos(sqrt((-4.*E*E*M+2*E*Q2+M*Q2)/(-4.*E*E*M+2*E*Q2)))
########################

########################
### Epsilon (polarization) for fixed value of Q2, initial energy E.
########################
def eps(E, Q2, M):
    return 1/(1 + 2.*(1+tau(Q2))/(4.*E*E/Q2-2.*E/M-1))
########################


########################
### 2. FORM FACTORS AND RADII
########################

########################
### Value of dipole form factor for fixed Q2, characteristic scale L2 = Lambda^2.
########################
def dipff(Q2, L2):
    return 1/(1+Q2/L2)**2
########################

########################
### Value of form factor using z expansion for fixed z and coefficient values {a_k}, so that
### G(z) = a0 + a_k*z**k for k > 0.
### par = [a_1, ..., a_kmax], so that a_0 = G(Q2=0) - a_k*z0**k, where z0 = z(Q2=0).
########################
def zff(z, a0, par):
    value = a0
    zn = z
    for i in range(len(par)):
        value += par[i]*zn
        zn = z*zn
    return value
########################

########################
### Radius in fm given coefficient values of form factor using z expansion with fixed t0, tcut.
### par = [a_1, ..., a_kmax], G0 = G(Q2=0).
########################
def zrad(t0, tcut, G0, par):
    z0 = (1-sqrt(1-t0/tcut))/(1+sqrt(1-t0/tcut))
    dGq20 = 0
    zn0 = 1
    for i in range(len(par)):
        dGq20 += (i+1)*par[i]*zn0
        zn0 = zn0*z0
    if dGq20<0:
        rad = sqrt(-6./tcut/G0*sqrt(1-t0/tcut)/((1+sqrt(1-t0/tcut))**2)*dGq20)*GeVfm
    else:
        rad = 0.0
    return rad
########################

########################
### Radius-SQ in fm given coefficient values of form factor using z expansion with fixed t0, tcut.
### par = [a_1, ..., a_kmax], G0 = G(Q2=0).
########################
def zrad2(t0, tcut, par):
    z0 = (1-sqrt(1-t0/tcut))/(1+sqrt(1-t0/tcut))
    dGq20 = 0
    zn0 = 1
    for i in range(len(par)):
        dGq20 += (i+1)*par[i]*zn0
        zn0 = zn0*z0
    
    rad2 = (-6./tcut*sqrt(1-t0/tcut)/((1+sqrt(1-t0/tcut))**2)*dGq20)*GeVfm**2 #divided by G0 outside
    
    return rad2
########################

########################
### Slope (a1) for form factor using z expansion with fixed t0, tcut,
### given radius in fm and higher-order coefficients, {G(Q2=0), rad, a_2, ..., a_kmax}, where
### G0 = G(Q2=0), par = [a_2, ..., a_kmax].
########################
def slope(t0, tcut, G0, rad, par):
    if rad > 0:
        a1 = -G0*tcut/6.*(1+sqrt(1-t0/tcut))**2/sqrt(1-t0/tcut)*(rad/GeVfm)**2
    else: # NEG ERROR
        a1 = G0*tcut/6.*(1+sqrt(1-t0/tcut))**2/sqrt(1-t0/tcut)*(rad/GeVfm)**2
    zn0 = z0
    for i in range(len(par)):
        a1 -= (i+2)*par[i]*zn0
        zn0 = zn0*z0
    return a1
########################

########################
### Value of polynomial form factor for fixed Q2, coefficient values, and G0 = G(Q2=0), so that
### G(Q2) = G0 + a_k*Q2**k for k > 0.
### G0 = G(Q2=0), par = [a_1, ..., a_kmax].
########################
def polyff(Q2, G0, par): 
    value = G0
    q = Q2
    for i in range(len(par)):
        value += par[i]*q
        q = q*Q2
    return value    
########################

########################
### Radius in fm given coefficient values of polynomial or inverse polynomial form factor.
### G0 = G(Q2=0), par = [a_1, ..., a_kmax].
########################
def polyrad(G0, par):
    return sqrt(abs(6./G0*par[0]))*GeVfm
########################


########################
### 3. FORM FACTOR RATIOS AND REDUCED CROSS SECTIONS
########################

########################
### Compute reduced cross section for fixed Q2 (with ta=Q2/4/M**2) and electron scattering angle th
### for form factors using the z expansion with given electric and magnetic coefficient values.
### GE = a0 + a_k*z**k, GM = b0 + b_k*z**k, k > 0.
### parge = [a_1, ..., a_kmax], pargm = [b_1, ..., b_kmax].
########################
def redcs(Q2, ta, th, z, a0, parge, b0, pargm):
    ##! If desired, insert code for (inverse) polynomial form factors here.
    
    # z expansion
    ge = a0
    gm = b0
    zn = z
    for i in range(len(parge)):
        ge += parge[i]*zn
        gm += pargm[i]*zn
        zn = z*zn
    ge2 = ge*ge
    gm2 = gm*gm

    # Compute reduced Rosenbluth cross section with model form factors.
    xsec = (ge2+ta*gm2)/(1+ta)+2*ta*gm2*tan(th/2)**2
    
    return xsec
########################


########################
## Compute FF based on z-expension
########################
def getff(Q2, z, a0, parg):
    #z expension
    g = a0
    zn = z
    for i in range(len(parg)):
        g += parg[i] *zn
        zn=z*zn
    return g
########################


########################
### Compute ratio (GE/GE0)/(GM/GM0) for fixed Q2 or z of electric and magnetic form factors GE0 and GM0,
### scaled by the G(Q2=0), with GE0 = GE(Q2=0), GM0 = GE(Q2=0).
### GE = a0 + a_k*z**k, GM = b0 + b_k*z**k, k > 0.
### parge = [a_1, ..., a_kmax], pargm = [b_1, ..., b_kmax].
########################
def ffratio(Q2, z, GE0, a0, parge, GM0, b0, pargm):
    ##! If desired, insert code for (inverse) polynomial form factors here.
    
    # z expansion
    ge = a0
    gm = b0
    zn = z
    for i in range(len(parge)):
        ge += parge[i]*zn
        gm += pargm[i]*zn
        zn = z*zn
    
    return GM0*ge/(GE0*gm)
########################


########################
### 4. SUM RULES FOR Z EXPANSION.
### Deprecated, now directly invert the matrix using numpy.linalg.inv.
########################

########################
### Compute matrix of sum rules as outlined in notes. These constrain the 4 highest-order coefficients 
### by enforcing the expected dipole-like high Q^2 behaviour of the form factor, 
### for fixed z0 and nmax = kmax - 4.
########################
def sumrules(z0, nmax):
    srmat = [[0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0], [0,0,0,0,0]]
    # Column 1
    srmat[0][0] = 6.
    srmat[1][0] = float(-(nmax+2)*(nmax+3)*(nmax+4))
    srmat[2][0] = float(3*(nmax+1)*(nmax+3)*(nmax+4))
    srmat[3][0] = float(-3*(nmax+1)*(nmax+2)*(nmax+4))
    srmat[4][0] = float((nmax+1)*(nmax+2)*(nmax+3))
    # Column 2
    srmat[0][1] = float(z0**(nmax+1)*((nmax+1)*(nmax+2)*(nmax+3)*z0**3 - 3*(nmax+1)*(nmax+2)*(nmax+4)*z0**2 + 3*((nmax+8)*nmax + 19)*nmax*z0 + 36*z0 - nmax*((nmax+9)*nmax + 26) - 24))
    srmat[1][1] = float((nmax+2)*(nmax+3)*(nmax+4))
    srmat[2][1] = float(-3*(nmax+1)*(nmax+3)*(nmax+4))
    srmat[3][1] = float(3*(nmax+1)*(nmax+2)*(nmax+4))
    srmat[4][1] = float(-(nmax+1)*(nmax+2)*(nmax+3))
    # Column 3
    srmat[0][2] = float(-3*(z0-1)*z0**(nmax+1)*(nmax**2 + (nmax+1)*(nmax+2)*z0**2 - 2*(nmax+1)*(nmax+3)*z0 + 5*nmax + 6))
    srmat[1][2] = float((nmax+2)*(nmax+3)*(z0**(nmax+2)*(nmax**2 + (nmax+1)*(nmax+2)*z0**2 - 2*(nmax+1)*(nmax+4)*z0 + 7*nmax + 12) - 6))/2
    srmat[2][2] = float(-(nmax+3)*(z0**(nmax+1)*(nmax*(nmax**2 + 2*(nmax+1)*(nmax+2)*z0**3 - 3*(nmax+1)*(nmax+4)*z0**2 + 9*nmax + 26) + 24) - 6*(3*nmax+4)))/2
    srmat[3][2] = float((nmax+1)*(z0**(nmax+1)*((nmax+2)*(nmax+1)*nmax*z0**3 - 3*nmax*(nmax+3)*(nmax+4)*z0 + 2*(nmax+2)*(nmax+3)*(nmax+4)) - 6*(3*nmax+8)))/2
    srmat[4][2] = float((nmax+1)*(nmax+2)*(6 - z0**(nmax+1)*(nmax*(z0-1)*(nmax*(z0-1) + z0 - 5) + 6)))/2
    # Column 4
    srmat[0][3] = float(3*(z0-1)**2*z0**(nmax+1)*(nmax*(z0-1) + z0 - 2))
    srmat[1][3] = float((nmax+2)*(3 - z0**(nmax+2)*(nmax**2*(z0-1)**2 + nmax*(4*z0-7)*(z0-1) + 3*(z0-2)**2)))
    srmat[2][3] = float(z0**(nmax+1)*(nmax**3 + 9*nmax**2 + (nmax+1)*(2*nmax+1)*(nmax+3)*z0**3 - 3*(nmax+1)**2*(nmax+4)*z0**2 + 26*nmax + 24) - 3*(3*nmax + 5))
    srmat[3][3] = float(-(2*nmax**3 + 15*nmax**2 + (nmax+2)*(nmax+1)*nmax*z0**3 - 3*(nmax+1)**2*(nmax+4)*z0 + 34*nmax + 24)*z0**(nmax+1) + 9*nmax + 12)
    srmat[4][3] = float((nmax+1)*(z0**(nmax+1)*(nmax*(z0-1)*(nmax*(z0-1) + 2*z0 - 5) - 3*z0 + 6) - 3))
    # Column 5
    srmat[0][4] = float(-(z0-1)**3*z0**(nmax+1))
    srmat[1][4] = float(z0**(nmax+2)*(nmax**2 + (nmax+2)*(nmax+3)*z0**2 - 2*(nmax+2)*(nmax+4)*z0 + 7*nmax + 12) - 2)/2
    srmat[2][4] = float(6 - z0**(nmax+1)*(nmax**2 + 2*(nmax+1)*(nmax+3)*z0**3 - 3*(nmax+1)*(nmax+4)*z0**2 + 7*nmax + 12))/2
    srmat[3][4] = float(z0**(nmax+1)*((nmax+1)*(nmax+2)*z0**3 - 3*(nmax+1)*(nmax+4)*z0 + 2*(nmax+2)*(nmax+4)))/2 - 3.
    srmat[4][4] = float(2 - z0**(nmax+1)*(nmax**2 + (nmax+1)*(nmax+2)*z0**2 - 2*(nmax+1)*(nmax+3)*z0 + 5*nmax + 6))/2
    # Normalize
    srmatnorm = -(nmax+2)*(nmax+3)*(nmax+4)*z0**(nmax+1) + 3*(nmax+1)*(nmax+3)*(nmax+4)*z0**(nmax+2) - 3*(nmax+1)*(nmax+2)*(nmax+4)*z0**(nmax+3) + (nmax+1)*(nmax+2)*(nmax+3)*z0**(nmax+4) + 6
    for i in range(len(srmat)):
        for j in range(len(srmat[i])):
            srmat[i][j] = srmat[i][j]/srmatnorm

    return srmat
########################
