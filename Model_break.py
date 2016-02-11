import sncosmo
import astropy
import os
import numpy as np
import matplotlib.pyplot as plt

'''
Note, the following are a small sample of the available built in models in sncosmo:
TypeIb model: source=sncosmo.get_source('snana-2006ep', version='1.0')
TypeIb/c model: sncosmo.get_source('nugent-sn1bc', version='1.1')
TypeIc model: sncosmo.get_source('snana-2004fe',version='1.0')
TypeIIP model: sncosmo.get_source('s11-2004hx',version='1.0')
'''

# Registering ptf48r filter
import os
import sys

sys.path.append(os.path.abspath('..'))
import Bandpass


sourcelist= [[ 'nugent-sn1bc', 'nugent-hyper','s11-2005hl', 's11-2005hm',  's11-2006jo','snana-2004gv', 'snana-2006ep', 'snana-2007y', 'snana-2004ib', 'snana-2005hm', 'snana-2007nc',],[ 'nugent-sn1bc', 'nugent-hyper', 's11-2006fo','snana-2004fe', 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo', 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la', 'snana-04d4jv',],[ 'nugent-sn2p', 's11-2004hx','s11-2005lc','s11-2005gi','s11-2006jl', 'snana-2007ms','snana-2004hx', 'snana-2005gi', 'snana-2006gq', 'snana-2006kn', 'snana-2006jl', 'snana-2006iw', 'snana-2006kv', 'snana-2006ns', 'snana-2007iz', 'snana-2007nr', 'snana-2007kw', 'snana-2007ky', 'snana-2007lj', 'snana-2007lb', 'snana-2007ll', 'snana-2007nw', 'snana-2007ld', 'snana-2007md', 'snana-2007lz', 'snana-2007lx', 'snana-2007og', 'snana-2007ny', 'snana-2007nv', 'snana-2007pg',],['nugent-sn2l', 'nugent-sn2n', 'snana-2006ez', 'snana-2006ix', 'whalen-z15b', 'whalen-z15d', 'whalen-z15g', 'whalen-z25b', 'whalen-z25d', 'whalen-z25g']]
#'snana-2006jo' 'whalen-z40b', 'whalen-z40g' removed


# t0 can either be the start time or the peak time. Each in built model has it's limits. The code determines
# at what times the model begins to break down, therefore determining whether t0 is peak time or start time.
def model_limits(sou):
    source = sncosmo.get_source(sou)
    model = sncosmo.Model(source=source)

    timearray = range(-20, 300)
    fluxlc = model.bandflux('ptf48r',timearray)     # Creating a flux array of the light curve
    modelbreak = []
    for i in range(len(fluxlc)-1):
        if fluxlc[i] == fluxlc[i+1]:
            modelbreak.append(timearray[i+1])
    for i in range(len(modelbreak)-1):
        if (modelbreak[i+1]-modelbreak[i]) > 20:
            return sou, modelbreak[i], modelbreak[i+1]


for sou in sourcelist[0][:]:
    try:
        print model_limits(sou)
    except:
        pass