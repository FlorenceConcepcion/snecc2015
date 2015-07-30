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
	source=sncosmo.get_source('snana-2004hx',version='1.0') 
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
	
        wave = np.arange(3000.,8000., 0.1)
        spec = model.flux([30.,40.,50.,60.], wave)
        print spec
        print maglc
        plt.plot(wave, spec[0], color = 'purple')
        plt.plot(wave, spec[1], color = 'blue')
        plt.plot(wave, spec[2], color = 'green')
        plt.plot(wave, spec[3], color = 'red')

                                                
#	plt.scatter(timearray, maglc2, color = 'red')

        #plt.ylim([0,3e-17])
	#plt.gca().invert_yaxis()
	plt.draw()
	
	'''
	m[:,2]	| ccdid
	m[:,3]	| lmt_mg_new
	m[:,4]	| Seeing_ratio
	m[:,5]	| medsky_new
	m[:,6]	| good_pix_area	
	'''
	'''
	##Getting Rid of NANs in the Mag arrays
	time_array=m[:,0][~np.isnan(maglc)]
	mag_lc=maglc[~np.isnan(maglc)]
	flux_lc=fluxlc[~np.isnan(maglc)]
	ccd_lc=m[:,2][~np.isnan(maglc)]
	lmt_lc=m[:,3][~np.isnan(maglc)]
	see_rat=m[:,4][~np.isnan(maglc)]
	med_lc=m[:,5][~np.isnan(maglc)]
	pix_lc=m[:,6][~np.isnan(maglc)]
	#print maglc
	#print mag_lc
	sn_par=np.array((time_array, mag_lc, flux_lc, ccd_lc, lmt_lc, see_rat, med_lc, pix_lc )).T
	'''
	'''snapr
	snpar[:,0]	| time
	snpar[:,1]	| mag_lc
	snpar[:,2]	| flux_lc
	snpar[:,3]	| ccd_lc
	snpar[:,4]	| lmt_lc
	snpar[:,5]	| see_rat
	snpar[:,6]	| med_lc
	snpar[:,7]	| pix_lc
	'''
	#print sn_par

	return  maglc, fluxlc, timearray
	
Gen_SN()
