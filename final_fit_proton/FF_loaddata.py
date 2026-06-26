######## 
### FF_loaddata.py
### Functions that extract relevant subsets of Mainz, world, and polarization datasets.
########
from math import *

########################
### CONSTANTS
########################
pi = 3.14159265358
alpha = 7.29735e-3 # Fine-structure constant.
Mp = 0.938272 # Mass of the proton in GeV.
mup = 2.792847356 # G_M(0), value of the Sachs magnetic form factor at q^2 = 0, which for the proton is the magnetic moment.
Lambda2 = 0.71 # Best fit value for parameter in dipole form factor for proton.
GeVfm = 0.197327 # GeV^(-1) to fm conversion: 1 GeV^(-1) = 0.197327 fm.
barn = 0.389379338e-3 # 1 GeV^(-2)=0.389e-3 barn.
########################


########################
### Load Mainz rebinned data with Q2 <= Q2max.
### Compute respective z(Q2; tcut, t0) for each Q2 point.
### A1opt is a flag to incorporate A1 systematics (different energy cuts or correlated systematics from kinematics).
### tpeopt is a flag to use either Feshbach or no TPE correction (instead of default Blunden).
########################
def loadMainz(Q2max, t0, tcut, A1opt, tpeopt):# {{{
    data = []
    normlist = [] # list of normalization combinations
    specElist = [[], [], []] # list of spectrometer-Ebeam combinations.
    ##! mainzscale[spec][Ebeam]: Mainz scaling factors for pt-pt errors, as applied by A1 collaboration detailed in Table I of Explanations.pdf.
    ##! mainzscale = [[1.070,1.753,1.654,1.215,1.465,1.100], [1.266,2.080,1.550,1.081,1.700,1.700], [1.366,2.293,1.833,1.766,1.100,1.170]]
    
    #for l in open('data/proton/Mainz_CrossSections_Rebinned.dat'):
    for l in open('data/proton/Mainz_CrossSections_Rebinned.dat'):
        values = l.split()
        
        # Primary kinematic variables are averaged Q^2 [GeV^2], beam energy [GeV].
        try: # eliminates descriptions in 1st line of file
            Q2 = float(values[3])
            en = float(values[0])/1000
        except:
            continue
        if en == 0.315 and float(values[2]) == 30.01: # remove point with large scatter
            continue
        
        # Q2max cut (or fixed energy).
        if Q2 > Q2max:
            continue

        # Spectrometer and normalization parameters.
        if values[1] == 'A': 
            spec = 0
        elif values[1] == 'B':
            spec = 1
        else:
            spec = 2
        norms = [int(i) for i in values[-1].split(":")] # normalization parameters
        if norms not in normlist:
            normlist.append(norms)
        if [spec, en] not in specElist[spec]:
            specElist[spec].append([spec, en])

        # Compute secondary kinematic quantities.
        ta = Q2/4./Mp**2 # tau = Q^2/(4*Mp**2)
        th = 2*acos(sqrt((-4.*en*en*Mp+2*en*Q2+Mp*Q2)/(-4.*en*en*Mp+2*en*Q2))) # electron scattering angle
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut)) # z
        Ef = en - Q2/2/Mp # final electron energy
            
        # Compute dipole reduced cross sections.
        ge2 = (1/(1+Q2/Lambda2)**2)**2 # square of dipole form factor.
        gm2 = ge2*mup**2
        dipcs = (ge2+ta*gm2)/(1+ta)+2*ta*gm2*tan(th/2)**2

        # TPE correction factor.
        if tpeopt == 1: # Feshbach
            tpe = float(values[-3])
        elif tpeopt == 2: # none
            tpe = 1
        else: # Blunden, default
            tpe = float(values[-4])
                
        # Expt reduced cs and uncertainty.
        dcs = float(values[5])/tpe*dipcs
        # Use different cs from A1 energy cut.
        if A1opt == 1: # smaller cut in energy loss, col 8
            cs = float(values[7])/tpe*dipcs
        elif A1opt == 2: # larger cut in energy loss, col 9
            cs = float(values[8])/tpe*dipcs
        else:
            cs = float(values[4])/tpe*dipcs

        # Apply correlated systematics from A1 analysis.
        if A1opt == 3: # Multiply cs by A1 systematic scale in col 12.
            cs = cs*float(values[-2])
            dcs = dcs*float(values[-2])
        elif A1opt == 4: # Divide cs by A1 systematic scale in col 12.
            cs = cs/float(values[-2])
            dcs = dcs/float(values[-2])

        # Append all the above as an array to data:
        # 0. energy, 1. Q^2, 2. tau, 3. theta, 4. z,
        # 5. expt. reduced xsec, 6. uncertainty, 7. normalization parameters, 8. spectrometer
        data.append([en, Q2, ta, th, z, cs, dcs, norms, spec])            

    return data, normlist, specElist
######################### }}}

########################
### Load world data with Q2 <= Q2max.
### Compute respective z(Q2; tcut, t0) for each Q2 point.
### tpeopt is a flag to use either Feshbach or no TPE correction (instead of default Blunden).
########################
def loadworld(Q2max, t0, tcut, tpeopt, ex_tpe_opt):# {{{
    data = []
    normlist = [] # list of experiment numbers and corresponding normalization uncertainties
    
    for l in open('data/proton/World_CrossSections_2017-03-30.dat'):
    #! for l in open('data/proton/World_CrossSections_Aug23.dat'):
        values = l.split()

        # Q2max cut.
        try: # eliminates descriptions in 1st line of file
            Q2 = float(values[2]) # Q2 [GeV^2] from file
        except:
            continue
        if Q2 > Q2max:
            continue

        # Primary kinematic variables are averaged energy [GeV] and electron scattering angle theta.
        en = float(values[0])/1000.
        th = float(values[1])/180.*pi

        # Experiment number and normalization uncertainty.
        expnum = int(values[7].rstrip('.')) # experiment number
        dnorm = float(values[8]) # experimental normalization uncertainty
        if [expnum, dnorm] not in normlist: # insert at beginning of normlist, as data file is sorted in descending order of exptnum 
            normlist.insert(0, [expnum, dnorm])

        # Compute other kinematic quantities.
        Q2 = 2*en**2/(1/(1-cos(th)) + en/Mp) # computed Q2 [GeV^2] instead of from file
        ta = Q2/4./Mp**2 # tau = Q^2/(4*Mp**2)
        #eps = 1/(1 + 2*(1+ta)/(4*en**2/Q2-2*en/Mp-1))
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut)) # z
        Ef = en/(1+en/Mp*(1-cos(th))) # final electron energy

        # Compute Mott cross section in nb/sr with recoil factor.
        mott = ((alpha*cos(th/2))/(2*en*sin(th/2)**2))**2*Ef/en*1e9*barn

        # TPE correction factor.
        if tpeopt == 1: # Feshbach
            tpe = float(values[6])
        elif tpeopt == 2: # none
            tpe = 1
        else: # Blunden, default
            tpe = float(values[5])

        # Expt reduced cs and uncertainty.
        cs = float(values[3])/tpe/mott
        dcs = float(values[4])/tpe/mott

        # Extra TPE at high Q2
        # ex_tpe_opt = 1
        if(ex_tpe_opt and Q2>1.0):
            eps = (1.0+2.0*(1.0+Q2/4.0/(Mp**2)) *(tan(th/2.))**2 )**(-1.)
            delta_extra = -0.01 * (1.0-eps) * log(Q2)/log(2.2)
            cs *= (1.0-delta_extra)
            dcs = cs* sqrt( (dcs/cs)**2 + (2.0*delta_extra)**2  )

        # Append all the above as an array to data:
        # 0. energy, 1. Q^2, 2. tau, 3. theta, 4. z,
        # 5. expt. reduced xsec, 6. uncertainty of (5.), 7. expt no, 8. expt norm uncertainty.
        data.append([en, Q2, ta, th, z, cs, dcs, expnum, dnorm])
    
    return data, normlist
########################
# }}}

#######################
### Load polarization data with Q2 <= Q2max.
### Compute respective z(Q2; tcut, t0) for each Q2 point.
########################
def loadpol(Q2max, t0, tcut):            # {{{
    data = []
    explist = [] # list of experiment numbers

    for l in open('data/proton/World_GeGm_Polarization.dat'):
    #! for l in open('data/proton/World_GEpGMp_Polarization.dat'):
        values = l.split()
            
        # Q2max cut.
        try: # eliminates descriptions in 1st line of file
            Q2 = float(values[0]) # Q2 [GeV^2]
        except:
            continue
        if Q2 > Q2max:
            continue

        # Experiment number.
        expnum = int(values[4])
        if expnum not in explist:
            explist.append(expnum)

        # Computed kinematic quantities.
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut)) # z

        # Ge/Gm and uncertainty.
        rat = float(values[1])
        drat = sqrt(float(values[2])**2+float(values[3])**2) # quadrature sum of stat and syst errors
        
        # Append as an array to datapol:
        # 0. Q^2, 1. z, 2. mup*GE/GM ratio, 3. uncertainty, 4. experiment number
        data.append([Q2, z, rat, drat, expnum])

    return data, explist
######################### }}}


########################
### Load fake data with Q2 <= Q2max., added by Z. Ye 08/26/2016
### Compute respective z(Q2; tcut, t0) for each Q2 point.
########################
def loadfakeHQ(Q2max, t0, tcut):            # {{{
    data = []
    # filename = 'data/proton/prot_GEGM_fakeHQ_sep20.dat'
    filename = 'data/proton/prot_GEGM_fakeHQ_sep20.dat'

    for l in open(filename):
        values = l.split()
            
        # Q2max cut.
        try: # eliminates descriptions in 1st line of file
            Q2 = float(values[0]) # Q2 [GeV^2]
        except:
            continue
        if Q2 > Q2max:
            continue

        # Computed kinematic quantities.
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut)) # z

        # Ge/Gm and uncertainty.
        GD = float(values[7])
        GE = float(values[1]) * GD
        dGE = float(values[2]) * GD
        GM = float(values[3]) * GD *mup
        dGM = float(values[4]) * GD *mup
        GEGM = float(values[5])
        dGEGM = float(values[6])
        
        # Append as an array to datafake
        data.append([Q2, z, GE, dGE, GM, dGM, GEGM, dGEGM])

    return data
######################### }}}

########################
### Load SBS polarization-ratio data with Q2 <= Q2max.
### Compute z(Q2; tcut, t0) for each point.
### Input file format (whitespace separated):
###   Q2   muGE/GM   d(muGE/GM)
########################
def loadSBS(Q2max, t0, tcut):                  # {{{
    data = []
    filename = 'data/GEpGMp_SBS.dat'

    for l in open(filename):
        l = l.strip()
        if (not l) or l.startswith('#'):
            continue

        values = l.split()
        if len(values) < 3:
            continue

        Q2    = float(values[0])
        GEGM  = float(values[1])
        dGEGM = float(values[2])

        if Q2 > Q2max:
            continue

        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut))

        # 0. Q^2, 1. z, 2. muGE/GM ratio, 3. uncertainty
        data.append([Q2, z, GEGM, dGEGM])

    return data
######################### }}}

########################
## Load GE&GM world data
## from John's paper Phy. Rev. C 76 035205(2007)
########################
def loadGEGM(Q2_max):# {{{
    Q2=[]
    GEdip=[]
    dGEdip=[]
    GMdip=[]
    dGMdip=[]
    GEGMdip=[]
    dGEGMdip=[]
    
    fitlines = open('data/proton/GEGM_John.dat', 'r').readlines()
    for i in range(1,len(fitlines)-1):
        values = fitlines[i].split()
        try:#eliminates descriptions in 1st line of file
            Q2test=float(values[0])
        except:
            continue

        if Q2test>Q2_max:
            continue

        Q2.append(values[0])
        GEdip.append(values[1])
        dGEdip.append(values[2])
        GMdip.append(values[3])
        dGMdip.append(values[4])
        GEGMdip.append(values[5])
        dGEGMdip.append(values[6])

    return Q2, GEdip, dGEdip, GMdip, dGMdip, GEGMdip, dGEGMdip
######################## }}}
