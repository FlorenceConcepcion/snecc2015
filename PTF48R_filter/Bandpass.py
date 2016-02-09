import sncosmo  
import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('../PTF48R_filter/PTF48R.dat')
wavelength = np.array([row[0] for row in data])
transmission = np.array([row[1] for row in data])
'''
plt.scatter(wavelength, transmission)
plt.show()
'''
band = sncosmo.Bandpass(wavelength, transmission, name='PTF48R')
sncosmo.registry.register(band)
#print 'transmission-weighted effective wavelength:', band.wave_eff
