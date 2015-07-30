#10qqd

import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt


def linefit(x):
    hml=np.load('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/' + x)    

    hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec'),
                                          dtype=('str','float','float','float','float','float','float','float','float','float'))
    
    x = np.array(hml[:,1], dtype = float)
    y = np.array(hml[:,2], dtype = float)
    yerr = np.array(hml[:,3], dtype = float)
    
    
    #averagey = np.mean(y)*np.ones(len(y))
    #print averagey
    plt.errorbar(x,y,yerr, fmt= 'o')
    plt.axhline(y = np.mean(y), color = 'purple')
    #plt.show()
    
    ChiSq = 0.0
    for i in hml[:,2]:
        ChiSq += ((float(i) - np.mean(y))**2/ np.var(y))
    
    #returns reduced Chi squared
    return ChiSq/float(len(hml[:,1])-1)


cclist = []
for i in os.walk('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                print linefit(k)
                cclist.append(k)
