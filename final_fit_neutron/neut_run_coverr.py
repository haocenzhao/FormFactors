#!/usr/bin/python3
# Filename: run_coverr.py
# This code runs neut_FF_sumrules_coverr.py
# and writes the results to output files.

import sys, os, traceback
import numpy as n

## General parameters.
bnd = '5'

kmin = 11
kmax = 11
Q2_plot_max = 1000.0

Gmode = int(sys.argv[1])
# max=4 for GEn, max=10.0 for GMn
gmod = 'GEn'
Q2min = 4
Q2max = 4
if Gmode == 2:
    gmod = 'GMn'
    Q2min = 10
    Q2max = 10

if Gmode == 1:
    gmod = 'GEn'
    Q2min = 10
    Q2max = 10

#t0mode = int(sys.argv[2])
t0mode = 4
if t0mode == 1:
    mod = 't0zero'
if t0mode == 2:
    mod = 't0opt'
if t0mode == 3:
    mod = 't0fix'
if t0mode == 4:
    mod = 't0fix7'
if t0mode == 5:
    mod = 't0fix4'

fitmode = int(sys.argv[2])
#fitmode = 1

Q2max = 1000
Q2min = Q2max
#Q2max = int(sys.argv[4])

for Q2 in range(Q2min, Q2max+1):
    for k in range(kmin, kmax+1):
        folder = 'Apr2026_world_norm_bound5/z' + str(k)

        fitdata = 'Apr2026_world'
        if fitmode == 1:
            input_file = folder + '/' + fitdata + '_sumrules_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_' + gmod + '.dat'
        if fitmode == 2:
            input_file = folder + '/' + fitdata + '_sumrules_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_' + gmod + '_noHQ.dat'
        if fitmode == 3:
            input_file = folder + '/' + fitdata + '_sumrules_leastsq_Q2' + str(Q2) + '_z' + str(k) + '_gb' + bnd + '_' + mod + '_' + gmod + '_noRad.dat'

        output_file = input_file.replace('/Apr2026_world_', '/out_world_')

        if 0:
            print('--- file has been created: %s, skipped!' % input_file)
            os.system('python3 neut_FF_sumrules_coverr.py %f %s %s' % (Q2_plot_max, input_file, output_file))
        elif os.path.isfile(input_file):
            print('--- Reading file = ', input_file)
            print('--- Output to file = ', output_file)
            os.system('python3 neut_FF_sumrules_coverr.py %f %s %s' % (Q2_plot_max, input_file, output_file))
        else:
            print('--- file can not be found: %s, skipped!' % input_file)

print('--- All done!!!')