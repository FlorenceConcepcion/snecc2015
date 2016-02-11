import sncosmo
import astropy
import os
import numpy as np
import matplotlib.pyplot as plt

# creating a dictionary with supernova name and type of supernova
x=open('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt')

'''
# Extracting data from a list of Supernovae, their name and type of sn
hml = open('List/SNlist.txt')
lines = [line.split(' ') for line in hml if len(line) > 2]
types = {}
for line in lines:
    for i in range(len(line)):
        if line[i] == 'SN':
            types[line[0]] = line[i + 1].strip(',')
'''
# Dictionary created from given list of core collapse supernova (using code above)
types = {'10hgi': 'Ib/c', '10pjg': 'IIP', '10qwz': 'IIP', '10rin': 'IIP', '10wmf': 'IIP', '10hdq': 'IIP',
         '10ttd': 'IIP', '10qqd': 'Ic', '10qts': 'Ic-BL', '10tpa': 'IIP', '10xfh': 'Ic', '10osr': 'IIP', '10pbi': 'Ib',
         '10wal': 'Ic', '10vgv': 'Ic-BL', '10hie': 'Ic', '10rmn': 'IIP', '10xik': 'Ic', '10nbf': 'IIP', '10glp': 'IIP',
         '10xem': 'Ic-BL', '10lbo': 'Ic', '10yow': 'Ic', '10osn': 'Ic', '10rem': 'IIP', '10hny': 'IIP', '10hpa': 'IIP',
         '10qif': 'Ib', '10tqv': 'Ic', '10vdl': 'IIP', '10uqn': 'IIP', '10inj': 'Ib', '10rjs': 'IIP', '10npd': 'IIP',
         '10hfe': 'Ic', '10qob': 'IIP', '10ood': 'Ic', '10uhf': 'Ic', '10svt': 'Ic', '10hyq': 'IIP', '10vwg': 'Ic',
         '10tqi': 'Ic', '10uhn': 'Ic', '10myz': 'IIP', '10wve': 'IIP', '10ysd': 'Ic-BL', '10vnv': 'Ib', '10hcw': 'Ib',
         '10tff': 'IIP', '10gtn': 'IIP', '10gvb': 'Ic-BL', '10kui': 'Ib', '10xjv': 'IIP', '10xjr': 'Ib',
         '10iue': 'Ic-BL'}

# List of models available in sncosmo
sourcelist= [[ 'nugent-sn1bc', 'nugent-hyper','s11-2005hl', 's11-2005hm',  's11-2006jo','snana-2004gv', 'snana-2006ep', 'snana-2007y', 'snana-2004ib', 'snana-2005hm', 'snana-2007nc',],[ 'nugent-sn1bc', 'nugent-hyper', 's11-2006fo','snana-2004fe', 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo', 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la', 'snana-04d4jv',],[ 'nugent-sn2p', 's11-2004hx','s11-2005lc','s11-2005gi','s11-2006jl', 'snana-2007ms','snana-2004hx', 'snana-2005gi', 'snana-2006gq', 'snana-2006kn', 'snana-2006jl', 'snana-2006iw', 'snana-2006kv', 'snana-2006ns', 'snana-2007iz', 'snana-2007nr', 'snana-2007kw', 'snana-2007ky', 'snana-2007lj', 'snana-2007lb', 'snana-2007ll', 'snana-2007nw', 'snana-2007ld', 'snana-2007md', 'snana-2007lz', 'snana-2007lx', 'snana-2007og', 'snana-2007ny', 'snana-2007nv', 'snana-2007pg',],['nugent-sn2l', 'nugent-sn2n', 'snana-2006ez', 'snana-2006ix', 'whalen-z15b', 'whalen-z15d', 'whalen-z15g', 'whalen-z25b', 'whalen-z25d', 'whalen-z25g']] 
#'snana-2006jo' 'whalen-z40b', 'whalen-z40g' removed


#this function requires the model name, z and t0 values so that the values returned are specific to each data set
def modellimit(source_name, z, t0):
    source=sncosmo.get_source(source_name) 
    model=sncosmo.Model(source=source)
    
    #Creating a flux array of the model lightcurve
    timearray= np.arange(t0-120,t0+500, 0.01)
    model.set(z=z,t0=t0)
    fluxlc=model.bandflux('ptf48r',timearray) 
    modelbreak = [t0-120.0]
    for i in range(len(fluxlc)-1):
        if fluxlc[i] == fluxlc[i+1]:
            modelbreak.append(timearray[i+1])

    #find the limits of the model by defining where the flux does not change. 
    for i in range(len(modelbreak)-1):
        if abs(modelbreak[i+1]-modelbreak[i]) > 20:
            minimum = t0 - modelbreak[i]
            maximum = modelbreak[i+1] - t0
            
            #return the difference of t0 and the edges of the model and the time of peak flux
            ##this is done as some models use t0 as peak and others use it for start of model
            return minimum, maximum, timearray[list(fluxlc).index(max(fluxlc))]

def fitcurve(x, source_name, hml):
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 
       # print len(hml[:,1]), 'data points'

        # Adding zpsys and filter columns. zp system used is ab and filter used is ptf48r
        ab = np.zeros(len(hml[:, 1]), dtype='|S2')
        band = np.zeros(len(hml[:, 1]), dtype='|S6')
        for i in range(len(ab)):
            ab[i] = 'ab'
            band[i] = 'ptf48r'
        hml = np.column_stack((hml, ab, band))

        # Converting into a table using astropy with titles: ptfname, time,
        # magnitude, mag_err, flux, flux_err, zp_new, zp, ra, dec, zpsys and filter
        hml_dat = astropy.table.Table(data=hml, names=(
             'ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys',
               'filter'),
                                      dtype=(
                                          'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float',
                                            'float', 'str', 'str'))
        # Fitting model: model parameter z bound.
        res, fitted_model = sncosmo.fit_lc(hml_dat, model, ['z', 't0', 'amplitude'], bounds={'z': (0.005, 0.35)})
        #The following excludes data points not in the range of the model and data sets with fewer than 4 data points
        limit = modellimit(source_name, res.parameters[0], res.parameters[1])
        hml2 = []
        for j in range(len(hml[:,1])):
            datapoint = hml [:,1][j]
            if (res.parameters[1]- limit[0])< float(datapoint) < (res.parameters[1]+limit[1]):
                hml2.append(hml[j]) 
        hml2 = np.array(hml2)
        if len(hml2)>3:
            return finalfitcurve(x, source_name, hml2)    
    
    except ValueError:
        pass   
         
def finalfitcurve(x, source_name, hml):
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 
        
        #fit to model
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(0.005,0.35)}, nburn=10000, nsamples=50000)
        
        #The following excludes data with not enough distribution..
        ##Ensuring at least 2 data points on either side of peak.
        limit = modellimit(source_name, res.parameters[0], res.parameters[1])
        j=0; k=0
        for i in hml[:,1]:
            if float(i)>float(limit[2]):
                j+=1
            else:
                k+=1        
        if j>=2 and k>=2:          
            return np.array([float(res.chisq/res.ndof), hml[:,0][0], source_name]), hml #returns reduced chi squared, sn name and model name. 
           
    except ValueError:
        pass

def showcurve(sn, source_name, hml):
    try:


        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 

        #fit to model
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(0.005,0.35)}, nburn=10000, nsamples=50000)

        #plot model
        sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name + '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=10)
        plt.show()

        #print 'Parameters:''z:',float(res.parameters[0]), float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])'
        print 'Done:', hml[:,0][0], 'z:',float(res.parameters[0]), 'Reduced chi^2:', res.chisq/res.ndof,#'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof

    except ValueError:
        print sn, source_name, 'cannot be plotted'

'''
# Lists all supernova data sets available in CorecollapseLC/ccdata/
cclist = []
for i in os.walk(datadir):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
print len(cclist), 'supernovae in list'
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
cclist = cclist_retrieved
            
testlist = cclist[10:]
testlist = ['10svt_LC.npy']
for sn in testlist:
    indexvalue = 0
    if types[sn[:-7]][1] == 'I':
        indexvalue = 2
    elif types[sn[:-7]][1] == 'b':
        indexvalue = 0
    else:
        indexvalue = 1
    fitted=[]; redchisq=[]
    #find sn data set to every model of sn type and record reduced chi squared
    for source in sourcelist[indexvalue][:1]:
        y = fitcurve(sn, source, np.load( datadir + sn))
        if y != None:
            fitted.append(y)
            redchisq.append(float(y[0][0]))
            
    if len(fitted)!=0:
        # Find model with best fit determined by minimised reduced chi squared
        minredchisq = min(redchisq)
        bestfit = fitted[redchisq.index(minredchisq)]
        
        #display plot
        showcurve(sn, bestfit[0][2],bestfit[1])

