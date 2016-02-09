# -*- coding: utf-8 -*-
import subprocess, math, sncosmo
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# Range of SN typeIa absolute magnitudes


alpha = 0.141;
beta = 3.101
values = []
for i in range(10000):
    x1 = np.random.uniform(-3, 3)
    c = np.random.uniform(-0.3, 0.3)
    sigma = np.random.normal(scale=0.15)
    absmag = -19.05 - alpha * x1 + beta * c + sigma
    values.append(absmag)

# plt.hist(values, bins = 100)
# plt.show()

# plotting a range of supernova by varying absolute magnitude

##Use below if not on iridis
source = sncosmo.get_source('salt2', version='2.4')
model = sncosmo.Model(source=source)
timearray = range(100)
mabs = -19.3  # Setting the absolute Magnitude of the Supernova
model.set(z=0.01, t0=30)  # Setting redshift

band = sncosmo.get_bandpass('ptf48r')  # Retrieving the ptf48r bandpass

for i in range(50):
    model.set_source_peakabsmag(values[i], 'ptf48r', 'ab')  # Fixing my peak absolute magnitude

    maglc = model.bandmag('ptf48r', 'ab', timearray)  # Creating a magnitude array of the lightcurve
    # maglc2=model.bandmag('sdssr','ab',timearray) #Creating a magnitude array of the lightcurve

    plt.scatter(timearray, maglc, color='blue', label='ptf48r')
    # plt.scatter(timearray, maglc2, color = 'red', label='sdssr')

# fluxlc=model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
# fluxlc2=model.bandflux('bessellr',timearray) #Creating a magnitude array of the lightcurve

plt.gca().invert_yaxis()
plt.title('SNTypeIa peak 30 days')
# plt.legend()
plt.show()
