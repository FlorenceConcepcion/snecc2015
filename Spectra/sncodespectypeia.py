# -*- coding: utf-8 -*-
import subprocess, math, sncosmo
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
plt.ion()
def Gen_SN():
	#Use below if on Iridis
	#source = sncosmo.SALT2Source(modeldir="/scratch/cf5g09/Monte_Carlos/salt2-4")

	##Use below if not on iridis
	source=sncosmo.get_source('hsiao',version='3.0') 
	model=sncosmo.Model(source=source) 
        timearray= range(100)
	mabs= -16.7  #Setting the absolute Magnitude of the Supernova
	model.set(z=0.01,t0=30) #Setting redshift
	
	model.set_source_peakabsmag(mabs,'bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
	#model.set(x1=x_1, c=colour)
	
	absmagb =model.source_peakabsmag('bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
	absmag_r =model.source_peakabsmag('bessellr','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))


	band=sncosmo.get_bandpass('bessellr') #Retrieving the ptf48r bandpass 
    

	maglc=model.bandmag('bessellb','ab',timearray) #Creating a magnitude array of the lightcurve  
	fluxlc=model.bandflux('bessellb',timearray) #Creating a flux array of the lightcurve
	maglc2=model.bandmag('bessellr','ab',timearray) #Creating a magnitude array of the lightcurve 
	
        wave = np.arange(3000.,9000., 0.1)
        spec = model.flux([20.,30.,35.,40.,50.], wave)

        plt.plot(wave, spec[1], color = '#2980b9',label = 'Peak')
        plt.plot(wave, spec[2], color = '#27ae60',label = '5 days after peak')
        plt.plot(wave, spec[3], color = '#c0392b',label = '10 days after peak')
        plt.plot(wave, spec[4], color = '#8e44ad',label = '20 days after peak')

        plt.bar(5277, max(spec[1])*1.2,width=300, color = '#9b59b6',alpha = .3, edgecolor = '#9b59b6')
        plt.bar(4000, max(spec[1])*1.2,width=80, color = '#2ecc71',alpha = .3,edgecolor='#2ecc71')
        plt.bar(6130, max(spec[1])*1.2,width=120, color = '#2ecc71',alpha = .3, edgecolor='#2ecc71')
        plt.bar(4700, max(spec[1])*1.2,width=400, color = '#3498db',alpha = .3,edgecolor = '#3498db')
        
        font=  {'family':'Lato', 'size': '15'}
        
        plt.text(3895, 2.4e-15, 'Silicon', fontweight = 'semibold')        
        plt.text(4822, 2.e-15, 'Iron', fontweight = 'semibold')
        plt.text(5263, 2.e-15, 'Sulphur', fontweight = 'semibold')
        plt.text(6045, 2.e-15, 'Silicon', fontweight = 'semibold')

        plt.title('Model spectrum of Type Ia supernova')
        plt.ylabel('Flux in ergs / s / cm^2 / $\AA$')
        plt.xlabel('Wavelength in $\AA$')
        plt.minorticks_on()
        plt.rc('font', **font)
        plt.ylim([0,max(spec[1])*1.2])
	plt.legend()
	plt.draw()

	return  maglc, fluxlc, timearray
	
Gen_SN()
