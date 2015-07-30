import subprocess, math, sncosmo, astropy
import numpy as np
import matplotlib.pyplot as plt

source=sncosmo.get_source('salt2',version='2.4') 
model=sncosmo.Model(source=source) 

try:
    hml=np.load('/home/fcm1g13/Documents/Supernova/Fit 10hml_LC/10hml_LC.npy')    
    print len(hml[:,1]), 'data points'
    
    #adding zpsys and filter columns
    ab=np.zeros(len(hml[:,1]), dtype='|S2')
    for i in range(len(ab)):
        ab[i]='ab'
    hml=np.column_stack((hml,ab))   
    band=np.zeros(len(hml[:,1]), dtype='|S6')
    for i in range(len(band)):
        band[i]='ptf48r'
    hml=np.column_stack((hml,band))
    hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str','str'))

    #fitting model
    print model
    res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','x0','x1','c'], bounds={'z':(0.005,0.35),'x1':(-3.5,3.5), 'c':(-0.35,0.35)}, nburn=100, nsamples=5000)
    print 'chi^2 value at minimum', res.chisq
    print ' number of chi^2 function calls made:', res.ncall
    print 'dof:', res.ndof
    print 'reduced chi^2:', res.chisq/res.ndof
    sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]), xfigsize=10)
    plt.show()
    print #'### Parameters ###'
    print 'SN',str(hml[:,0][0]), 'z:',float(res.parameters[0])# float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])
    print 'Done:', hml[:,0][0]
        
except ValueError:
    print 'Value Error' 
        
        