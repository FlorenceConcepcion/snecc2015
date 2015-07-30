
import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM


phase = np.linspace(-50., 50., 5)
disp = np.linspace(5000., 9000., 6)
x = np.array([[0.], [2.], [4.], [2.], [0.]])
flux = np.repeat(x, 6, axis=1)

flux[0]=flux[0]*10
flux[1]=flux[1]*10

for z in flux:
    z[0]+=10
    z[1]+=5

print phase
print disp
print flux
#flux = np.transpose(flux)
model = sncosmo.TimeSeriesSource(phase, disp, flux)

timearray= np.arange(-120,120, 0.01)
#mabs= -16.7  #Setting the absolute Magnitude of the Supernova
#model.set(z=0.01,t0=30) #Setting redshift

#model.set_source_peakabsmag(mabs,'bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
#model.set(x1=x_1, c=colour)
#band=sncosmo.get_bandpass('bessellr') #Retrieving the bandpass 
#maglc=model.bandmag('bessellb','ab',timearray) #Creating a magnitude array of the lightcurve  
fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
fluxlc2=model.bandflux('bessellr',timearray) #Creating a magnitude array of the lightcurve 
fluxlc3=model.bandflux('sdssr',timearray) #Creating a magnitude array of the lightcurve 

'''
modelbreak = [-1000.0]
for i in range(len(fluxlc)-1):
    if fluxlc[i] == fluxlc[i+1]:
        modelbreak.append(timearray[i+1])
for i in range(len(modelbreak)-1):
    if (modelbreak[i+1]-modelbreak[i]) > 20:
        print modelbreak[i], modelbreak[i+1]

'''

plt.scatter(timearray, fluxlc, color= 'blue')
plt.scatter(timearray, fluxlc2, color = 'red')
plt.scatter(timearray, fluxlc3, color = 'green')
'''
spec = model.flux(phase, disp)

for i in range(len(phase)):
    print phase
    print spec[i]
    plt.plot(disp, spec[i], label = phase[i])

#plt.gca().invert_yaxis()
plt.legend()
'''
plt.show()





