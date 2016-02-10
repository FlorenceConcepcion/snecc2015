import astropy, sncosmo
import matplotlib.pyplot as plt
import numpy as np

'''
This code fits a Type Ia supernova model(10hml) to a Type Ia sn data set (salt2).
This code is simple and assumes there is enough data per data set to make a 
sensible fit and has bound values of z etc. which must be collected from
other sources. This code has limited variables and is built for a 'complete'
data set (as this supernova, 10hml, was well documented).
'''

# Registering ptf48r filter
import os
import sys
sys.path.append(os.path.abspath('..'))
import Bandpass

# Retrieving model from source, (requires internet access)
source = sncosmo.get_source('salt2', version='2.4')  # model used is salt2
model = sncosmo.Model(source=source)

# Load file containing data from same folder where code is stored
hml = np.load('10hml_LC.npy')
# print len(hml[:,1]), 'data points'         #number of rows and therefore data points

# Adding zpsys and filter columns. zp system used is ab and filter used is ptf48r
ab = np.zeros(len(hml[:, 1]), dtype='|S2')
band = np.zeros(len(hml[:, 1]), dtype='|S6')
for i in range(len(ab)):
    ab[i] = 'ab'
    band[i] = 'ptf48r'
hml = np.column_stack((hml, ab, band))

# Converting into a table using astropy with titles: ptfname, time,
# magnitude, mag_err, flux, flux_err, zp_new, zp, ra, dec, zpsys and filter
hml_dat = astropy.table.Table(data=hml, names=(
    'ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'),
                              dtype=(
                                  'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',
                                  'float', 'str', 'str'))

# Fitting model: model parameters of z, x1 and c are bound.
res, fitted_model = sncosmo.fit_lc(hml_dat, model, ['z', 't0', 'x0', 'x1', 'c'],
                                   bounds={'z': (0.005, 0.35), 'x1': (-3.5, 3.5), 'c': (-0.35, 0.35)})

# Use sncosmo to plot data and error
sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:, 0][0]), xfigsize=10)

# The following information can be shown if wished  :
# print 'chi^2 value at minimum:', res.chisq, ';', 'dof:', res.ndof
print 'Number of chi^2 function calls made:', res.ncall
print 'reduced chi^2:', res.chisq / res.ndof
print  # '### Parameters ###'
print 'SN', str(hml[:, 0][0]), '; z:', float(res.parameters[0])
print 'Completed fit for supernova:', hml[:, 0][0]

# Displays plot
plt.show()