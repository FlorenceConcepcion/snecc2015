# -*- coding: utf-8 -*-
import subprocess, math, sncosmo
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
plt.ion()
def Gen_SN():

	##Use below if not on iridis
	source=sncosmo.get_source('snana-2006gq')
	model=sncosmo.Model(source=source) 
        timearray= np.arange(-100,90, 0.01)
	mabs= -16.7  #Setting the absolute Magnitude of the Supernova
	model.set(z=0.01,t0=0) #Setting redshift
	
	#model.set_source_peakabsmag(mabs,'bessellb', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
	#model.set(x1=x_1, c=colour)
	
	#absmagb =model.source_peakabsmag('bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
	#absmag_r =model.source_peakabsmag('bessellr','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))


	#band=sncosmo.get_bandpass('bessellr') #Retrieving the ptf48r bandpass 
    

	maglc=model.bandmag('bessellb','ab',timearray) #Creating a magnitude array of the lightcurve  
	fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
	fluxlc2=model.bandflux('sdssg',timearray) #Creating a magnitude array of the lightcurve 
	fluxlc3=model.bandflux('sdssr',timearray) #Creating a magnitude array of the lightcurve 
	fluxlc4=model.bandflux('sdssi',timearray) #Creating a magnitude array of the lightcurve 
	fluxlc5=model.bandflux('sdssz',timearray) #Creating a magnitude array of the lightcurve 
	#fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
	
	modelbreak = [0.0]
	for i in range(len(fluxlc)-1):
	    if fluxlc[i] == fluxlc[i+1]:
	        modelbreak.append(timearray[i+1])
	for i in range(len(modelbreak)-1):
	    if (modelbreak[i+1]-modelbreak[i]) > 20:
	        print modelbreak[i], modelbreak[i+1]

        
	plt.plot(timearray, fluxlc, color= '#8e44ad', label = 'ptf48')
	plt.plot(timearray, fluxlc2, color = '#27ae60', label = 'Green (g)')
	plt.plot(timearray, fluxlc3, color = '#e74c3c', label = 'Red (r)')
	plt.plot(timearray, fluxlc4, color = '#34495e', label = 'Near Infrared (i)')
	plt.plot(timearray, fluxlc5, color = 'black', label = 'Infrared (z)')


	#plt.gca().invert_yaxis()

	plt.title('Light curve of type Ia model as seen in various sdss filters')
	plt.xlabel('Days after peak')
	plt.ylabel('Relative')
	plt.legend()
	#plt.yscale('log')
	plt.show()
	

	#print sn_par

	#return  maglc, fluxlc, timearray
	
Gen_SN()
