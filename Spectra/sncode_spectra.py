import sncosmo
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

'''
The aim of this code is to display the observed spectrum of a supernova at time intervals and to
identify peaks/troughs caused by specific elements. To see the spectra for each supernova type,
please change the source and comment out the code for other types of supernovae. This are also
saved as .png files.

Note, the following are a small sample of the available built in models in sncosmo:
TypeIa model: source=sncosmo.get_source('hsiao',version='3.0')
TypeIb model: sncosmo.get_source('snana-2007y',version='1.0')
TypeIc model: sncosmo.get_source('snana-2004fe',version='1.0')
TypeIIP model: sncosmo.get_source('snana-2004hx',version='1.0')
'''


def Gen_SN():
    # Retrieving model from source, (requires internet access)
    source = sncosmo.get_source('snana-2004hx',version='1.0')
    model = sncosmo.Model(source=source)

    # Setting absolute Magnitude, redshift and time of the Supernova.
    mabs = -16.7
    model.set(z=0.01, t0=30)  # (t0 is either start time or peak time depending on model.)
    model.set_source_peakabsmag(mabs, 'bessellr', 'ab', cosmo=FlatLambdaCDM(H0=70, Om0=0.3))

    # Plotting relative flux of different times against flux
    wave = np.arange(3000., 9000., 0.1)
    spec = model.flux([20., 30., 35., 40., 50.], wave)
    font = {'family': 'Lato', 'size': '15'}

    '''
    The following is supernova type specific, and to some extent specific to the model used.
    Please Note: To plot them all as subplots, the source and model must be changed for every plot.
    '''
    '''
    #-----------------------------------------------------------------
    # Model spectrum of Type Ia supernova
    plt.title('Model spectrum of Type Ia supernova')
    plt.plot(wave, spec[1], color='#2980b9', label='Peak')
    plt.plot(wave, spec[2], color='#27ae60', label='5 days after peak')
    plt.plot(wave, spec[3], color='#c0392b', label='10 days after peak')
    plt.plot(wave, spec[4], color='#8e44ad', label='20 days after peak')

    plt.bar(5277, max(spec[1]) * 1.2, width=300, color='#9b59b6', alpha=.3, edgecolor='#9b59b6')
    plt.bar(4000, max(spec[1]) * 1.2, width=80, color='#2ecc71', alpha=.3, edgecolor='#2ecc71')
    plt.bar(6130, max(spec[1]) * 1.2, width=120, color='#2ecc71', alpha=.3, edgecolor='#2ecc71')
    plt.bar(4700, max(spec[1]) * 1.2, width=400, color='#3498db', alpha=.3, edgecolor='#3498db')

    plt.text(3895, 2.4e-15, 'Silicon', fontweight='semibold')
    plt.text(4822, 2.e-15, 'Iron', fontweight='semibold')
    plt.text(5263, 2.e-15, 'Sulphur', fontweight='semibold')
    plt.text(6045, 2.e-15, 'Silicon', fontweight='semibold')

    #-----------------------------------------------------------------
    # Model spectrum of Type Ib supernova
    plt.title('Model spectrum of Type Ib supernova')
    plt.plot(wave, spec[1], color = '#2980b9',label = 'Peak')
    plt.plot(wave, spec[2], color = '#27ae60',label = '5 days after peak')
    plt.plot(wave, spec[3], color = '#c0392b',label = '10 days after peak')
    plt.plot(wave, spec[4], color = '#8e44ad',label = '20 days after peak')

    plt.bar(6230, max(spec[1])*1.2,width=100, color = '#c0392b',alpha = .3, edgecolor = '#c0392b')
    plt.bar(4800, max(spec[1])*1.2,width=130, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
    plt.bar(5682, max(spec[1])*1.2,width=90, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
    plt.bar(6494, max(spec[1])*1.2,width=60, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
    plt.bar(6860, max(spec[1])*1.2,width=100, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
    plt.bar(4050, max(spec[1])*1.2,width=500, color = '#3498db',alpha = .3, edgecolor = '#3498db')

    plt.text(6000, 2.03e-15, 'Hydrogen', fontweight = 'semibold')
    plt.text(4745, 1.8e-15, 'Helium', fontweight = 'semibold')
    plt.text(5580, 1.8e-15, 'Helium', fontweight = 'semibold')
    plt.text(6400, 1.8e-15, 'Helium', fontweight = 'semibold')
    plt.text(6800, 1.8e-15, 'Helium', fontweight = 'semibold')
    plt.text(4220, 2.03e-15, 'Iron', fontweight = 'semibold')

    #-----------------------------------------------------------------
    # Model spectrum of Type Ic supernova
    plt.title('Model spectrum of Type Ic supernova')
    plt.plot(wave, spec[1], color = '#2980b9',label = 'Peak')
    plt.plot(wave, spec[2], color = '#27ae60',label = '5 days after peak')
    plt.plot(wave, spec[3], color = '#c0392b',label = '10 days after peak')
    plt.plot(wave, spec[4], color = '#8e44ad',label = '20 days after peak')

    plt.bar(3980, max(spec[1])*1.2,width=560, color = '#3498db',alpha = .3, edgecolor = '#3498db')
    plt.bar(4643, max(spec[1])*1.2,width=590, color = '#3498db',alpha = .3,edgecolor = '#3498db')
    plt.bar(6125, max(spec[1])*1.2,width=150, color = '#2ecc71',alpha = .3,edgecolor = '#2ecc71')

    plt.text(4150, 2.03e-15, 'Iron', fontweight = 'semibold')
    plt.text(4833, 2.03e-15, 'Iron', fontweight = 'semibold')
    plt.text(6070, 2.03e-15, 'Silicon', fontweight = 'semibold')
    '''
    #-----------------------------------------------------------------
    # Model spectrum of Type IIP supernova
    plt.title('Model spectrum of Type IIP supernova')
    plt.plot(wave, spec[1], color = '#2980b9',label = 'Peak')
    plt.plot(wave, spec[2], color = '#27ae60',label = '5 days after peak')
    plt.plot(wave, spec[3], color = '#c0392b',label = '10 days after peak')
    plt.plot(wave, spec[4], color = '#8e44ad',label = '20 days after peak')

    plt.bar(4240, max(spec[1])*1.1,width=60, color = '#3498db',alpha = .3, edgecolor = '#3498db')
    plt.bar(6550, max(spec[1])*1.1,width=100, color = '#e74c3c',alpha = .3, edgecolor = '#e74c3c')
    plt.bar(4750, max(spec[1])*1.1,width=50, color = '#e74c3c',alpha = .3,edgecolor = '#e74c3c')

    plt.text(4200, 2.5e-15, 'Iron', fontweight = 'semibold')
    plt.text(6420, 2.5e-15, 'Hydrogen', fontweight = 'semibold')
    plt.text(4575, 2.5e-15, 'Hydrogen', fontweight = 'semibold')

    #-----------------------------------------------------------------
    plt.ylabel('Flux in ergs / s / cm^2 / $\AA$')
    plt.xlabel('Wavelength in $\AA$')
    plt.minorticks_on()
    plt.rc('font', **font)
    plt.ylim([0, max(spec[1]) * 1.2])

    plt.legend()
    plt.show()

Gen_SN()


