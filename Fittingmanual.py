
import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

cclist = []
for i in os.walk('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
#print cclist



plt.ion()
phase = np.linspace(2455454.88645-50., 2455454.88645+50., 5)
disp = np.linspace(1000., 9000., 6)
x = np.array([[0.], [2.], [4.], [2.], [0.]])
flux = np.repeat(x, 6, axis=1)

flux[0]=flux[0]*10
flux[1]=flux[1]*10

for z in flux:
    z[0]+=10
    z[1]+=5
    z[4]+=5
    z[5]+=10    

print phase
print disp
print flux
flux = flux*0.0000000000002
#flux = np.transpose(flux)
model = sncosmo.TimeSeriesSource(phase, disp, flux)

timearray= np.arange(2455454.88645-120., 2455454.88645+120, 0.01)
#mabs= -16.7  #Setting the absolute Magnitude of the Supernova
#model.set(z=0.01,t0=30) #Setting redshift

#model.set_source_peakabsmag(mabs,'bessellb','ab', cosmo=FlatLambdaCDM(H0=70,Om0=0.3)) #Fixing my peak absolute magnitude
#model.set(x1=x_1, c=colour)
#band=sncosmo.get_bandpass('bessellr') #Retrieving the bandpass 
#maglc=model.bandmag('bessellb','ab',timearray) #Creating a magnitude array of the lightcurve  
fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
fluxlc2=model.bandflux('bessellr',timearray) #Creating a magnitude array of the lightcurve 
fluxlc3=model.bandflux('sdssr',timearray) #Creating a magnitude array of the lightcurve 

plt.scatter(timearray, fluxlc, color= 'blue')
#plt.scatter(timearray, fluxlc2, color = 'red')
#plt.scatter(timearray, fluxlc3, color = 'green')


hml = np.load('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/10osr_LC.npy')
#testlist = ['10svt_LC.npy']

ab=np.zeros(len(hml[:,1]), dtype='|S2')
for i in range(len(ab)):
    ab[i]='ab'
hml=np.column_stack((hml,ab))   
band=np.zeros(len(hml[:,1]), dtype='|S6')
for i in range(len(band)):
    band[i]='ptf48r'
hml=np.column_stack((hml,band))
#fit to model
z0 = 0.001
hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))


res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['amplitude'], bounds={}, nburn=10000, nsamples=50000)

#plot model

#sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name + '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=10)
#plt.show()

#print 'Parameters:''z:',float(res.parameters[0]), float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])'
print 'Done:', hml[:,0][0], 'z:',float(res.parameters[0]), 'Reduced chi^2:', res.chisq/res.ndof,#'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof

#plt.errorbar(np.array(hml[:,1], dtype=float), np.array(hml[:,4], dtype=float), yerr = np.array(hml[:,5], dtype=float), fmt='o', color = 'red')

#plt.show()
