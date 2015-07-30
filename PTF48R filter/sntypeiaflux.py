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
	source=sncosmo.get_source('salt2',version='2.4') 
	model=sncosmo.Model(source=source) 
        timearray= range(100)
	mabs= -19.05  #Setting the absolute Magnitude of the Supernova
	model.set(z=0.01,t0=30) #Setting redshift
	
	model.set_source_peakabsmag(mabs,'ptf48r','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
	#model.set(x1=x_1, c=colour)
	
	absmag_r =model.source_peakabsmag('ptf48r','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
        print absmag_r

	band=sncosmo.get_bandpass('ptf48r') #Retrieving the ptf48r bandpass 
    

	maglc=model.bandmag('ptf48r','ab',timearray) #Creating a magnitude array of the lightcurve  
	fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve


	model2 = model
	model2.set_source_peakabsmag(mabs,'bessellr','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
	#model.set(x1=x_1, c=colour)
	
	absmag_r =model2.source_peakabsmag('bessellr','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3))
	maglc2=model2.bandmag('bessellr','ab',timearray) #Creating a magnitude array of the lightcurve  

	fluxlc2=model2.bandflux('bessellr',timearray) #Creating a magnitude array of the lightcurve 
	


        plt.scatter(timearray, fluxlc, color= 'blue', label='ptf48r')
	plt.scatter(timearray, fluxlc2, color = 'red', label='bessellr')

	model.bandflux('PTF48R', 30)
	


	#plt.gca().invert_yaxis()
	plt.title('SNTypeIa peak 30 days')
	plt.legend()
	plt.show()
	

	#print sn_par

	return  maglc, fluxlc, timearray
	
Gen_SN()
