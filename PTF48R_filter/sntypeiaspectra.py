# -*- coding: utf-8 -*-
import math
import sncosmo
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM

plt.ion()


def Gen_SN():
    source = sncosmo.get_source('salt2', version='2.4')
    model = sncosmo.Model(source=source)
    timearray = range(100)
    mabs = -19.05  # Setting the absolute Magnitude of the Supernova
    model.set(z=0.01, t0=30)  # Setting redshift

    model.set_source_peakabsmag(mabs, 'bessellb', 'ab',
                                cosmo=FlatLambdaCDM(H0=70, Om0=0.3))  # Fixing my peak absolute magnitude
    # model.set(x1=x_1, c=colour)
    '''
    absmagb =model.source_peakabsmag('bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
    absmag_r =model.source_peakabsmag('bessellr','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))


    band=sncosmo.get_bandpass('bessellr') #Retrieving the ptf48r bandpass
    

    maglc=model.bandmag('bessellb','ab',timearray) #Creating a magnitude array of the lightcurve
    fluxlc=model.bandflux('bessellb',timearray) #Creating a flux array of the lightcurve
    maglc2=model.bandmag('bessellr','ab',timearray) #Creating a magnitude array of the lightcurve
    '''

    data = np.loadtxt('/home/fcm1g13/Documents/Supernova/PTF48R filter/PTF48R.dat')
    wavelength = np.array([row[0] for row in data])
    transmission = np.array([row[1] for row in data])

    band = sncosmo.Bandpass(wavelength, transmission, name='PTF48R')

    wave = np.arange(5580., 7500., 20)
    spec = model.flux([20., 30., 35., 40., 50.], wave)

    adjspec = transmission * spec[1]
    '''
        for point in transmission:
            for flux in spec[1]:
                adjspec.append(point*flux)
        '''
    print len(transmission)
    print len(spec[1])
    print adjspec

    plt.plot(wave, spec[1], color='#27ae60', label='Flux')
    model.bandflux('PTF48R', 30)
    plt.plot(wave, np.array(adjspec), color='#2980b9', label='Observed flux')
    plt.plot(wavelength, transmission * max(spec[1]))
    # plt.plot(wave, spec[2], color = '#27ae60',label = '5 days after peak')
    # plt.plot(wave, spec[3], color = '#c0392b',label = '10 days after peak')
    # plt.plot(wave, spec[4], color = '#8e44ad',label = '20 days after peak')

    plt.title('Model spectrum of Type Ia supernova at peak')
    plt.ylabel('Flux in ergs / s / cm^2 / $\AA$')
    plt.xlabel('Wavelength in $\AA$')
    plt.minorticks_on()
    plt.ylim([0, max(spec[1]) * 1.2])
    plt.legend()
    plt.draw()

    return timearray


Gen_SN()
