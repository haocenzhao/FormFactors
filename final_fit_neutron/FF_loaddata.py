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
### Load GEn/GMn data with Q2 <= Q2max.
### Compute respective z(Q2; tcut, t0) for each Q2 point.
########################
def loadGN(Q2max, t0, tcut, opt):            # {{{
    data = []
    explist = [] # list of experiment numbers
    normlist = [] # list of experiment numbers and corresponding normalization uncertainties
    input_data = ''
    if opt==1:
        input_data = './data/neutron/GLOBFIT17_gen_feb07.out'
    elif opt==2:
        input_data = './data/neutron/GLOBFIT17_gmn_feb07.out' 
    else:
        print("*** ERROR, I don't know the option(not 1 or 2) = ",opt)
        return data, explist
    print('-- Loading neutron FF data from %s'%input_data)

    for l in open(input_data):
        values = l.split()
            
        # Q2max cut.
        try: # eliminates descriptions in 1st line of file
            Q2 = float(values[0]) # Q2 [GeV^2]
        except:
            continue
        if Q2 > Q2max:
            continue

        # Experiment number.
        expnum = float(values[4])
        if opt==1:
            dnorm = 0.0 ## fix at 2% for GEn only
        if opt==2:
            dnorm = float(values[5]) # experimental normalization uncertainty
            
        if [expnum, dnorm] not in normlist: # insert at beginning of normlist, as data file is sorted in descending order of exptnum 
            if dnorm>1e-33:#remove 0 normalization value
                normlist.insert(0, [expnum, dnorm])
 
        if expnum not in explist:
            explist.append(expnum)

        # Computed kinematic quantities.
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut)) # z

        # Ge/Gm and uncertainty.
        rat = float(values[1])
        drat = float(values[2])
        
        # Append as an array to dataGN:
        # 0. Q^2, 1. z, 2. GN/GD ratio, 3. uncertainty, 4. experiment number
        #data.append([Q2, z, rat, drat, expnum])
        # 0. Q^2, 1. z, 2. GN/GD ratio, 3. uncertainty, 4. experiment number, 5. normalization
        data.append([Q2, z, rat, drat, expnum, dnorm])

    #return data, explist
    return data, explist, normlist
######################### }}}

########################
### Load fake data to constraint high Q2, added by Z. Ye 10/04/2016
### Compute respective z(Q2; tcut, t0) for each Q2 point.
########################
def loadgnfakeHQ(Q2max, t0, tcut, nopt):            # {{{
    data = []
    if nopt==1:
        filename = 'data/neutron/neut_GEn_fakeHQ.dat'; ##GEn
    if nopt==2:
        filename = 'data/neutron/neut_GMn_fakeHQ.dat'; ##GMn/GD/mu
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
        GN = float(values[1])  ##GEn or GMn/mu/GD
        dGN = float(values[2]) ##dGEn or dGMn/mu/GD
        
        # Append as an array to datafake
        data.append([Q2, z, GN, dGN])

    return data
######################### }}}

########################
### Load SBS GEn data with Q2 <= Q2max.
### Compute z(Q2; tcut, t0) for each point.
### Convert GE/GD to GE using GD from the input file.
########################
def load_SBS_GEn(Q2max, t0, tcut):            # {{{
    data = []
    filename = 'data/GEnGMn_SBS.dat'

    for l in open(filename):
        values = l.split()

        # Skip empty or comment lines.
        if len(values) == 0:
            continue
        if values[0].startswith('#'):
            continue

        # Expected columns:
        # 0: Q2, 1: GD, ..., 6: GEGD, 7: GEGD_Err
        try:
            Q2 = float(values[0])
            GD = float(values[1])
            GEGD = float(values[6])
            dGEGD = float(values[7])
        except:
            continue

        # Q2max cut.
        if Q2 > Q2max:
            continue

        # Compute z.
        z = (sqrt(1+Q2/tcut)-sqrt(1-t0/tcut))/(sqrt(1+Q2/tcut)+sqrt(1-t0/tcut))

        # Convert GE/GD to GE.
        GE = GEGD * GD
        dGE = dGEGD * GD

        # Append as an array to data:
        # 0. Q^2, 1. z, 2. GE, 3. uncertainty
        data.append([Q2, z, GE, dGE])

    return data
######################### }}}
