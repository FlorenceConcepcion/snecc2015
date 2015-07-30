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
	source=sncosmo.get_source('snana-2007y',version='1.0') 
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

        plt.bar(6230, max(spec[1])*1.2,width=100, color = '#c0392b',alpha = .3, edgecolor = '#c0392b')
        plt.bar(4800, max(spec[1])*1.2,width=130, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
        plt.bar(5682, max(spec[1])*1.2,width=90, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
        plt.bar(6494, max(spec[1])*1.2,width=60, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
        plt.bar(6860, max(spec[1])*1.2,width=100, color = '#f1c40f',alpha = .3, edgecolor = '#f1c40f')
        plt.bar(4050, max(spec[1])*1.2,width=500, color = '#3498db',alpha = .3, edgecolor = '#3498db')
 
        font=  {'family':'Lato', 'size': '15'}           

        plt.text(6000, 2.03e-15, 'Hydrogen', fontweight = 'semibold')        
        plt.text(4745, 1.8e-15, 'Helium', fontweight = 'semibold')
        plt.text(5580, 1.8e-15, 'Helium', fontweight = 'semibold')
        plt.text(6400, 1.8e-15, 'Helium', fontweight = 'semibold')
        plt.text(6800, 1.8e-15, 'Helium', fontweight = 'semibold')
        plt.text(4220, 2.03e-15, 'Iron', fontweight = 'semibold')
 
        plt.title('Model spectrum of Type Ib supernova')
        plt.ylabel('Flux in ergs / s / cm^2 / $\AA$')
        plt.xlabel('Wavelength in $\AA$')
        plt.minorticks_on()
        plt.rc('font', **font)
        plt.ylim([0,max(spec[1])*1.2])

	plt.legend()
	plt.draw()
                               
        #plt.ylim([0,3e-17])
	#plt.gca().invert_yaxis()
	plt.legend()
	plt.draw()

	#print sn_par

	return  maglc, fluxlc, timearray
	
Gen_SN()
