############################/
############################/
## Parameterized Form Factor Central Value and Error
############################/
## ID = 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn,
## Q2 in GeV^2
##
# The parameterization formula returns the uncertainty devided by G(0)*GD, where
#  GD(Q2) = 1./(1+Q2/0.71)^2
# and GEp(0) = 1, GMp(0) = 2.79284356, GEn(0) = 1, GMn(0) = -1.91304272,
#
# The parameterization formula for the Form Factor value is:
#  $$ GN(z) = sum_{i=0}^{N=12}(a_i * z^i)
# Note that the return value has been divided by (G(Q2=0)*G_Dip)
#
# The parameterization formula for the Form Factor error is:
# $$ log_{10}\frac{\delta G}{G_D} = (L+c_0)\Theta_a(L_1-L)
#                                 +\sum_{i=1}^{N}(c_i+d_i L)[\Theta_a(L_i-L)-\Theta_a(L_{i+1}-L)]
#                                 +log_{10}(E_{\inf})\Theta_a(L-L_{N+1})$$
# where $L=log_{10}(Q^2)$, $\Theta_{a}(x)=[1+10^{-ax}]^{-1}$. $a=1$.
import numpy as np
import copy
from math import *
from scipy.linalg import det

#################################################
#### Calculate the uncertainty for GMn from 
#### covariance matrix, accounting for sum rules
#################################################{{{

#############
# Solving for values with sum rules can easily just use the np.linalg package
# However, since we explicitly pass the Jacobian rather than having the fitter estimate it, we need to have an explicit form of the calculations
# Doing the calculations in the Jacobian function every time is quite computationally intensive
# Doing it by hand leaves _a lot_ of room for error
# This does some precalculations for Cramer's Rule of the pieces that are independent of the fit values so that number of operations are minimized
# Not as useful in this script, as it's not vectorized, but this function is useful for other scripts that are
#############
def sum_rule_determinant_precalc(N, tcut, t0):
    z0 = 0
    if t0 != 0:
        z0 = (sqrt(tcut+0.)-sqrt(tcut-t0))/(sqrt(tcut+0.)+sqrt(tcut-t0))

    nSumRules = 5
    sum_rule_matrix = np.array(
        [
            [1, z0 ** (N + 1), z0 ** (N + 2), z0 ** (N + 3), z0 ** (N + 4)],
            [1, 1, 1, 1, 1],
            [0, N + 1, N + 2, N + 3, N + 4],
            [0, N * (N + 1), (N + 1) * (N + 2), (N + 2) * (N + 3), (N + 3) * (N + 4)],
            [
                0,
                (N - 1) * N * (N + 1),
                N * (N + 1) * (N + 2),
                (N + 1) * (N + 2) * (N + 3),
                (N + 2) * (N + 3) * (N + 4),
            ],
        ]
    )

    # Determinant of sum rule matrix for denominator in Cramer's Rule
    detSRM = det(sum_rule_matrix)

    #############
    # Calculating each term requires substituting the column for the sum rule values
    # The Jacobian then requires we take the partial derivative of these calculations
    # Since the partial derivative is then dependent on a_k, we must have an explicit form
    # The sum rule values are the only piece here that is dependent on a_k
    # So, we can pre-calculate the coefficients of these and then simply multiply them by the partial derivative of the sum rule values in the Jacobian function
    #############

    # Holds the sum rule partial derivative coefficients
    sr_dets = np.zeros(shape=(nSumRules, nSumRules))

    # Loop over columns, corresponding to each variable
    for col in range(len(sum_rule_matrix)):
        # Don't change the original matrix
        num_mat1 = copy.deepcopy(sum_rule_matrix)

        # Zero out the column
        for row in range(len(sum_rule_matrix)):
            num_mat1[row][col] = 0

        # Loop over rows, corresponding to each sum rule
        for row in range(len(sum_rule_matrix)):
            # For easily unzeroing the row
            num_mat2 = copy.deepcopy(num_mat1)

            # Zero out the row
            for c2 in range(len(sum_rule_matrix)):
                num_mat2[row][c2] = 0

            # Only keep the pieces that will be multiplied by the sum rule partial derivative
            num_mat2[row][col] = 1

            sr_dets[row][col] = det(num_mat2) / detSRM

    return sr_dets

def poly_err_with_sum_rules(x, p, cov, tcut, t0):
    c = sum_rule_determinant_precalc(len(p), tcut=tcut, t0=t0)

    z0 = 0
    if t0 != 0:
        z0 = (sqrt(tcut+0.)-sqrt(tcut-t0))/(sqrt(tcut+0.)+sqrt(tcut-t0))

    N = len(p)

    jacobian = np.array([[(x**k) for k in range(1, N + 1)]])

    db_dak = [
        lambda k: z0**k,
        lambda k: 1,
        lambda k: k,
        lambda k: k * (k - 1),
        lambda k: k * (k - 1) * (k - 2),
    ]

    for i in range(5):
        for ji, j in enumerate([0, N + 1, N + 2, N + 3, N + 4]):
            jacobian = (
                jacobian
                + np.array(
                    [c[i, ji] * db_dak[i](k) * (x**j) for k in range(1, N + 1)]
                ).T
            )

    return np.sqrt(np.einsum("ij,jk,ki->i", jacobian, cov, jacobian.T))

def GetFF(kID, kQ2):# {{{
    ### GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
    if kID<1 or kID>4:
        # print '*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]'
        return -1000, -1000

    #################################################
    #### z-Expansion Parameters for Form Factor Values
    #################################################{{{
    GN_Coef_Fit = np.zeros((4,13), dtype=float)
    # GEp (SBS-projected update) 04/14/2026 by H. Zhao
    GN_Coef_Fit[0] = np.array([
        0.23964565944563468, -1.1121184847673562, 1.4483887337371693,
        0.5060816911674428, -2.3812423710630273, 0.964813042807631,
        2.2107798731715644, -4.474784240413616, 2.5123997623410066,
        4.523756974651614, -8.692343462041645, 5.488239443868501,
        -1.2336166229049361
    ])  # GEp

    # GMp (SBS-projected update) 04/14/2026 by H. Zhao
    MU_P = 2.792847356
    GN_Coef_Fit[1] = np.array([
        0.7388915687606409, -3.0660167270305467, 3.4308442078146686,
        1.7571197379823975, -3.8417489674962613, -3.3005271564231533,
        2.411294015856522, 13.179276106213102, -12.983897883078239,
        -13.389405287929634, 29.20924976651264, -18.104189475975772,
        3.9591100947936013
    ]) / MU_P  # GMp
    
    # GEn (SBS-projected update, z7) 04/14/2026 by H. Zhao 
    GN_Coef_Fit[2] = np.array([
        0.048838508462326874, -0.06491124817401979, -0.23631795283535903,
        0.4004663462195233, 0.25751343801333026, -0.6779328488081449,
        -0.015273646986079728, 0.5244143327404727, -0.24364568092946293,
        0.027674860377757682, -0.045679052826409716, 0.02485294474608013,
        0.0
    ]) #GEn
    
    #GMn -- REPLACED WITH NEW GLOBAL FIT - Oct 25 2025 ~TJH
    GN_Coef_Fit[3] = np.array([
        0.26056254, -1.09079786, 1.10759603, 
        1.13262363, -1.44919755, -3.9158823, 
        5.844183  ,  3.45835262, -11.79774606, 
        8.54072074, -2.09041479, 0.0, 0.0
    ]) #GMn
    #}}}
    
    GMN_fit_pars = np.load('fit_parSBS.npy')
    GMN_covmat = np.load('covmatSBS.npy')

#################################################
#### Parameters for Form Factor Errors
#################################################{{{
    parL = np.zeros((4,2), dtype=float)
    parM = np.zeros((4,15), dtype=float)
    parH = np.zeros((4,3), dtype=float)
    ## GEp:#!wrong
    parL[0] = np.array([-0.97775297,  0.99685273]) #Low-Q2
    parM[0] = np.array([ -1.97750308e+00,  -4.46566998e-01,   2.94508717e-01,   1.54467525e+00,
        9.05268347e-01,  -6.00008111e-01,  -1.10732394e+00,  -9.85982716e-02,
        4.63035988e-01,   1.37729116e-01,  -7.82991627e-02,  -3.63056932e-02,
        2.64219326e-03,   3.13261383e-03,   3.89593858e-04 ]) #Mid-Q2:
    parH[0] = np.array([ 0.78584754,  1.89052183, -0.4104746]) #High-Q2

    #GMp:!wrong
    parL[1] = np.array([-0.68452707,  0.99709151]) #Low-Q2
    parM[1] = np.array([ -1.76549673e+00,   1.67218457e-01,  -1.20542733e+00,  -4.72244127e-01,
        1.41548871e+00,   6.61320779e-01,  -8.16422909e-01,  -3.73804477e-01,
        2.62223992e-01,   1.28886639e-01,  -3.90901510e-02,  -2.44995181e-02,
        8.34270064e-04,   1.88226433e-03,   2.43073327e-04]) #Mid-Q2:
    parH[1] = np.array([  0.80374002,  1.98005828, -0.69700928]) #High-Q2
    
    #GEn:#!wrong
    parL[2] = np.array([-2.02311829, 1.00066282]) #Low-Q2
    parM[2] = np.array([-2.07343771e+00,   1.13218347e+00,   1.03946682e+00,  -2.79708561e-01,
        -3.39166129e-01,   1.98498974e-01,  -1.45403679e-01,  -1.21705930e-01,
        1.14234312e-01,   5.69989513e-02,  -2.33664051e-02,  -1.35740738e-02,
        7.84044667e-04,   1.19890550e-03,   1.55012141e-04,]) #Mid-Q2:
    parH[2] = np.array([0.4553596, 1.95063341, 0.32421279]) #High-Q2:

    GMnCovMat = np.array([[ 6.11351195e-05, -2.53001900e-05, -1.80527663e-04, -1.09073160e-04,  6.21116600e-04, -3.73268020e-04],
                 [-2.53001900e-05,  1.32692176e-04, -9.82614957e-06, -5.31334887e-04,  6.68574252e-04, -2.32305687e-04],
                 [-1.80527663e-04, -9.82614958e-06,  1.35868266e-03, -6.12499885e-04, -2.52584085e-03,  2.02553050e-03],
                 [-1.09073160e-04, -5.31334887e-04, -6.12499885e-04,  5.26981808e-03, -5.50107376e-03,  1.43297525e-03],
                 [ 6.21116600e-04,  6.68574252e-04, -2.52584085e-03, -5.50107376e-03,  1.36975403e-02, -7.04009888e-03],
                 [-3.73268020e-04, -2.32305687e-04,  2.02553050e-03,  1.43297525e-03, -7.04009888e-03,  4.26963335e-03]])

    ## Apply the z-expansion formula
    tcut = 0.0779191396
    t0 = -0.7
    z = (sqrt(tcut+kQ2)-sqrt(tcut-t0))/(sqrt(tcut+kQ2)+sqrt(tcut-t0)) 
    GNQ2 = 0
    if kID!=4:
        GNQ2 = np.array([GN_Coef_Fit[kID-1][i]*(z**i) for i in range(0, len(GN_Coef_Fit[kID-1]))]).sum() 
    else:
        GNQ2 = np.polynomial.polynomial.polyval(z, GMN_fit_pars).sum() 

    GDip= 1./(1. + kQ2/0.71)**2
    GNGD_Fit = GNQ2 / GDip #Note that the GN_Coef_Fit has been divided by mu_p or mu_n for GMp and GMn

    ## Apply the parameterization formula for error
    lnGNGD_Err=0.0
    GNGD_Err=0
    if kID!=4:
        lnQ2 = log10(kQ2)
        if kQ2<1e-3:
            lnGNGD_Err = parL[kID-1][0] + parL[kID-1][1]*lnQ2
        elif kQ2>1e2:
            lnGNGD_Err = parH[kID-1][0]*np.sqrt(lnQ2 - parH[kID-1][1]) + parH[kID-1][2]
        else:
            lnGNGD_Err = np.array([parM[kID-1][i]*(lnQ2**i) for i in range(0, len(parM[kID-1]))]).sum() 
        GNGD_Err = 10.**(lnGNGD_Err)    ##LOG10(dG/G(0)/GD)
    else:
        p_core = GMN_fit_pars[1:-4]
        cov_core = GMN_covmat[:-7, :-7]
        n = len(p_core)
        cov_core = cov_core[:n, :n]  # enforce shape (n, n)
        GNGD_Err = poly_err_with_sum_rules(z, p_core, cov_core, tcut, t0)[0] / GDip

    return GNGD_Fit, GNGD_Err
# }}}
