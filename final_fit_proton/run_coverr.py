#!/usr/bin/python3
# Filename: run_coverr.py
# This code runs FF_sumrules_coverr.py
# and writes the results to output files.

import sys, os, traceback
import numpy as n

## General parameters.
bnd = '5'
Q2_plot_max = 1000

Q2 = 1000
#Q2 = 1
#Q2 = 0.5

kmin = 12
kmax = 12

fittype = int(sys.argv[1])
outtype = int(sys.argv[2])
rad = 8409
radM = 851

fitdata = 'all'
if fittype == 1:
    fitdata = 'all'
if fittype == 2:
    fitdata = 'world'
if fittype == 3:
    fitdata = 'Mainz'

mod = 't0fix7'
#mod = 't0opt'

com = 'Apr2026'

for k in range(kmin, kmax+1):
    folder = '%s_%s_FF/z%s' % (com, fitdata, str(k))   # folder containing central fits
    fitdata1 = '%s_sumrules' % fitdata

    if outtype == 1:
        input_file = folder + '/' + fitdata1 + '_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '.dat'
    if outtype == 2:
        input_file = folder + '/' + fitdata1 + '_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_RE%d_RM%d.dat' % (rad, radM)
    if outtype == 3:
        input_file = folder + '/' + fitdata1 + '_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_noHQ.dat'
    if outtype == 4:
        input_file = folder + '/' + fitdata1 + '_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_noRad.dat'
    if outtype == 5:
        input_file = folder + '/' + fitdata1 + '_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_noExtraTPE.dat'

    output_file = input_file.replace('/%s_' % fitdata1, '/out_%s_' % fitdata1)
    print(' Input = ', input_file)
    print('Output = ', output_file)
    os.system('python3 FF_sumrules_coverr.py %f %s %s %d' % (Q2_plot_max, input_file, output_file, fittype))

print('--- All done!!!')