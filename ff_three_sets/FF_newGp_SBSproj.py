############################/
############################/
## Parameterized Form Factor Central Value and Error for Fake SBS Proton Data
## Edit by Haocen Zhao   01/27/2026
############################/
## ID = 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn,
## Q2 in GeV^2
##
# The parameterization returns:
#   - central value: G/(G(0)*GD)
#   - uncertainty:   dG/(G(0)*GD)
# where GD(Q2) = 1./(1+Q2/0.71)^2.
#
# Update in this file:
#   - Replace GEp (kID=1) and GMp (kID=2) z-expansion coefficients with SBS-projected fit.
#   - Keep GEn/GMn coefficients and all uncertainty parameterizations unchanged.

import numpy as np
from math import sqrt, log10


def GetFF(kID, kQ2):  # {{{
    """Return (G/(G0*GD), dG/(G0*GD)) for selected form factor ID at given Q2."""
    if kID < 1 or kID > 4:
        print('*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]')
        return -1000, -1000

    if kQ2 <= 0.0:
        print('*** ERROR***, Q2 must be positive.')
        return -1000, -1000

    #################################################
    # z-Expansion Parameters for Form Factor Values
    #################################################{{{
    GN_Coef_Fit = np.zeros((4, 13), dtype=float)

    # GEp (SBS-projected update)
    GN_Coef_Fit[0] = np.array([
        0.2394827349710868,  -1.1140157017539436,  1.4555556909100709,
        0.5389719001712749,  -2.474931783158097,   0.7886000931968082,
        2.74076680187188,    -4.226353939257555,   1.1538690881778528,
        5.0139701695575045,  -7.513591110692591,   4.3102937469430245,
        -0.912617690937287
    ])  # GEp

    # GMp (SBS-projected update)
    MU_P = 2.792847356
    GN_Coef_Fit[1] = np.array([
        0.7391443668755441,  -3.063388839383516,   3.424332588033044,
        1.7219230128791787,  -3.7678419618022447,  -3.1284540867860455,
        1.9972488932295311,  12.901803850365118,   -11.872206063178474,
        -13.693195077246315, 28.151823488855598,   -17.09892396807186,
        3.687733796230418
    ]) / MU_P  # GMp

    # GEn
    GN_Coef_Fit[2] = np.array([
        0.048919981379, -0.064525053912, -0.240825897382, 0.392108744873,
        0.300445258602, -0.661888687179, -0.175639769687, 0.624691724461,
        -0.077684299367, -0.236003975259, 0.090401973470, 0.0, 0.0
    ])  # GEn

    # GMn
    GN_Coef_Fit[3] = np.array([
        0.257758326959, -1.079540642058, 1.182183812195, 0.711015085833,
        -1.348080936796, -1.662444025208, 2.624354426029, 1.751234494568,
        -4.922300878888, 3.197892727312, -0.712072389946, 0.0, 0.0
    ])  # GMn
    #}}}

    #################################################
    # Parameters for Form Factor Errors 
    #################################################{{{
    parL = np.zeros((4, 2), dtype=float)
    parM = np.zeros((4, 15), dtype=float)
    parH = np.zeros((4, 3), dtype=float)

    # GEp
    parL[0] = np.array([-0.97775297, 0.99685273])#Low-Q2
    parM[0] = np.array([
        -1.97750308e+00, -4.46566998e-01, 2.94508717e-01, 1.54467525e+00,
        9.05268347e-01, -6.00008111e-01, -1.10732394e+00, -9.85982716e-02,
        4.63035988e-01, 1.37729116e-01, -7.82991627e-02, -3.63056932e-02,
        2.64219326e-03, 3.13261383e-03, 3.89593858e-04
    ])#Mid-Q2:
    parH[0] = np.array([0.78584754, 1.89052183, -0.4104746])#High-Q2

    # GMp
    parL[1] = np.array([-0.68452707, 0.99709151])#Low-Q2
    parM[1] = np.array([
        -1.76549673e+00, 1.67218457e-01, -1.20542733e+00, -4.72244127e-01,
        1.41548871e+00, 6.61320779e-01, -8.16422909e-01, -3.73804477e-01,
        2.62223992e-01, 1.28886639e-01, -3.90901510e-02, -2.44995181e-02,
        8.34270064e-04, 1.88226433e-03, 2.43073327e-04
    ])#Mid-Q2:
    parH[1] = np.array([0.80374002, 1.98005828, -0.69700928])#High-Q2

    # GEn
    parL[2] = np.array([-2.02311829, 1.00066282])
    parM[2] = np.array([
        -2.07343771e+00, 1.13218347e+00, 1.03946682e+00, -2.79708561e-01,
        -3.39166129e-01, 1.98498974e-01, -1.45403679e-01, -1.21705930e-01,
        1.14234312e-01, 5.69989513e-02, -2.33664051e-02, -1.35740738e-02,
        7.84044667e-04, 1.19890550e-03, 1.55012141e-04
    ])
    parH[2] = np.array([0.4553596, 1.95063341, 0.32421279])

    # GMn
    parL[3] = np.array([-0.20765505, 0.99767103])
    parM[3] = np.array([
        -2.07087611e+00, 4.32385770e-02, -3.28705077e-01, 5.08142662e-01,
        1.89103676e+00, 1.36784324e-01, -1.47078994e+00, -3.54336795e-01,
        4.98368396e-01, 1.77178596e-01, -7.34859451e-02, -3.72184066e-02,
        1.97024963e-03, 2.88676628e-03, 3.57964735e-04
    ])
    parH[3] = np.array([0.50859057, 1.96863291, 0.2321395])
    #}}}

    # z mapping
    tcut = 0.0779191396
    t0 = -0.7
    z = (sqrt(tcut + kQ2) - sqrt(tcut - t0)) / (sqrt(tcut + kQ2) + sqrt(tcut - t0))

    # Polynomial in z
    GNQ2 = np.array(
        [GN_Coef_Fit[kID - 1][i] * (z ** i) for i in range(len(GN_Coef_Fit[kID - 1]))]
    ).sum()

    # Divide by dipole
    GDip = 1.0 / (1.0 + kQ2 / 0.71) ** 2
    GNGD_Fit = GNQ2 / GDip

    # Error model in log10(Q2)
    lnQ2 = log10(kQ2)
    if kQ2 < 1e-3:
        lnGNGD_Err = parL[kID - 1][0] + parL[kID - 1][1] * lnQ2
    elif kQ2 > 1e2:
        lnGNGD_Err = parH[kID - 1][0] * np.sqrt(lnQ2 - parH[kID - 1][1]) + parH[kID - 1][2]
    else:
        lnGNGD_Err = np.array(
            [parM[kID - 1][i] * (lnQ2 ** i) for i in range(len(parM[kID - 1]))]
        ).sum()

    GNGD_Err = 10.0 ** lnGNGD_Err  ##LOG10(dG/G(0)/GD)

    return GNGD_Fit, GNGD_Err
# }}}