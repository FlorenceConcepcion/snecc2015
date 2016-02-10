# coding=utf-8
import subprocess, math, sncosmo
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

'''
Note, the following is a small sample of the available built in models in sncosmo:
salt2  - TypeIa model
s11-2005hm - TypeIb model
snana-2004fe TypeIc model
snana-2004hx - TypeIIP model
'''

def Gen_SN():
    # Retrieving model from source, (requires internet access)
    source = sncosmo.get_source('snana-2004hx', version='1.0')
    model = sncosmo.Model(source=source)


    # Registering ptf48r filter
    import os
    import sys
    sys.path.append(os.path.abspath('..'))
    import Bandpass

    # Setting absolute Magnitude, redshift and time of the Supernova.
    mabs = -16.7
    model.set(z=0.01, t0=30)  # (t0 is either start time or peak time depending on model.)
    model.set_source_peakabsmag(mabs, 'ptf48r', 'ab', cosmo=FlatLambdaCDM(H0=70, Om0=0.3))

    '''
    # As stated above t0 can either be the start time or the peak time. Each in built model has it's limits. The code
    # below can be used to determine at what times the model begins to break down, therefore determining whether
    # t0 is peak time or start time. 
    timearray = range(10, 100)
    fluxlc = model.bandflux('ptf48r',timearray) #Creating a flux array of the lightcurve
    modelbreak = [0.0]
    for i in range(len(fluxlc)-1):
        if fluxlc[i] == fluxlc[i+1]:
            modelbreak.append(timearray[i+1])
    for i in range(len(modelbreak)-1):
        if (modelbreak[i+1]-modelbreak[i]) > 20:
            print modelbreak[i], modelbreak[i+1]
    '''

    # Plotting relative flux of different wavelengths against time
    time_intervals = [25., 30., 40., 50.]
    wave = np.arange(3000., 8000., .1)
    spec = model.flux(time_intervals, wave)
    colors = ['black','red','green','blue']

    plt.subplot(122)
    for i in range(len(time_intervals)):
        plt.plot(wave, spec[i], color = colors[i%4], label =str(str(time_intervals[i]) + ' days'))
    plt.xlabel('Wavelength')
    plt.ylabel('Relative Intensity')
    plt.legend()

    # Plotting relative flux over time in different filters
    plt.subplot(121)
    timearray = range(10, 100)
    filters = ['ptf48r', 'bessellr', 'bessellv', 'bessellb' ]
    for i in range(len(time_intervals)):
        plt.scatter(timearray, model.bandflux(filters[i], timearray), color = colors[i%4], label= str(filters[i]))
    plt.xlabel('Time/ days')
    plt.ylabel('Relative Flux')
    plt.legend()

    plt.show()

    '''
    Magnitude can also be displayed using:
    plt.scatter(timearray, model.bandmag(filters[i], 'ab', timearray), color = colors[i%4], label= str(filters[i]))

    But, note, as it is magnitude, and so a logarithmic scale, there is less to note on a graph.
    '''

    # A 3D plot with axes representing time, wavelength and observed luminosity
    time_intervals = np.arange(23., 45., 1.)
    wave = np.arange(3000., 8000., 80.)
    spec = model.flux(time_intervals, wave)


    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    coords = []
    x = []
    y = []
    z = []
    for i in range(len(spec)):
        x.append([])
        y.append([])
        z.append([])
        for j in range(len(spec[0])):
            x[i].append(time_intervals[i])
            y[i].append(wave[j])
            z[i].append(spec[i][j])
        ax.scatter(x[i], y[i], z[i], c = [np.arange(0,0.9,0.2) for i in range(len(spec[0]))], marker='.')#cmap='nipy_spectral',

    ax.set_xlabel('Time /days')
    ax.set_ylabel('Wavelength /A')
    ax.set_zlabel('Relative Intensity')
    plt.show()


Gen_SN()

'''
# Lists all supernova data sets available in CorecollapseLC/ccdata/
cclist = []
for i in os.walk('CorecollapseLC/ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
print cclist
'''
cclist_retrieved = ['10glp_LC.npy', '10gtn_LC.npy', '10gvb_LC.npy', '10hcw_LC.npy', '10hdq_LC.npy', '10hfe_LC.npy',
                    '10hgi_LC.npy', '10hie_LC.npy', '10hny_LC.npy', '10hpa_LC.npy', '10hyq_LC.npy', '10inj_LC.npy',
                    '10iue_LC.npy', '10kui_LC.npy', '10lbo_LC.npy', '10myz_LC.npy', '10nbf_LC.npy', '10npd_LC.npy',
                    '10ood_LC.npy', '10osn_LC.npy', '10osr_LC.npy', '10pbi_LC.npy', '10pjg_LC.npy', '10qif_LC.npy',
                    '10qob_LC.npy', '10qqd_LC.npy', '10qts_LC.npy', '10qwz_LC.npy', '10rem_LC.npy', '10rin_LC.npy',
                    '10rjs_LC.npy', '10rmn_LC.npy', '10svt_LC.npy', '10tff_LC.npy', '10tpa_LC.npy', '10tqi_LC.npy',
                    '10tqv_LC.npy', '10ttd_LC.npy', '10uhf_LC.npy', '10uhn_LC.npy', '10uqn_LC.npy', '10vdl_LC.npy',
                    '10vgv_LC.npy', '10vnv_LC.npy', '10vwg_LC.npy', '10wal_LC.npy', '10wmf_LC.npy', '10wve_LC.npy',
                    '10xem_LC.npy', '10xfh_LC.npy', '10xik_LC.npy', '10xjr_LC.npy', '10xjv_LC.npy', '10yow_LC.npy',
                    '10ysd_LC.npy']

##But test the path out again