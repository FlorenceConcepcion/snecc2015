import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
import psycopg2



#creating a dictionary with supernova name and type of supernova  
x=open('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt')  
#fil = np.genfromtxt('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt', delimiter = '  ')#, dtype=(str, float, float,float, float,float, float,float,str,str,str,str))

lines = [line.split(' ') for line in x if len(line)>2]
types = {}

for line in lines:
    for i in range(len(line)):
        if line[i]=='SN':
            types[line[0]] = line[i+1].strip(',')
#list of models
sourcelist= [[ 'nugent-sn1bc', 'nugent-hyper','s11-2005hl', 's11-2005hm',  's11-2006jo','snana-2004gv', 'snana-2006ep', 'snana-2007y', 'snana-2004ib', 'snana-2005hm', 'snana-2007nc',],[ 'nugent-sn1bc', 'nugent-hyper', 's11-2006fo','snana-2004fe', 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo', 'snana-sdss014475', 'snana-2006lc', 'snana-04d1la', 'snana-04d4jv',],[ 'nugent-sn2p', 's11-2004hx','s11-2005lc','s11-2005gi','s11-2006jl', 'snana-2007ms','snana-2004hx', 'snana-2005gi', 'snana-2006gq', 'snana-2006kn', 'snana-2006jl', 'snana-2006iw', 'snana-2006kv', 'snana-2006ns', 'snana-2007iz', 'snana-2007nr', 'snana-2007kw', 'snana-2007ky', 'snana-2007lj', 'snana-2007lb', 'snana-2007ll', 'snana-2007nw', 'snana-2007ld', 'snana-2007md', 'snana-2007lz', 'snana-2007lx', 'snana-2007og', 'snana-2007ny', 'snana-2007nv', 'snana-2007pg',],['nugent-sn2l', 'nugent-sn2n', 'snana-2006ez', 'snana-2006ix', 'whalen-z15b', 'whalen-z15d', 'whalen-z15g', 'whalen-z25b', 'whalen-z25d', 'whalen-z25g']] 
#'snana-2006jo' 'whalen-z40b', 'whalen-z40g' removed


#get spectroscopic redshift from database
def spec_z(sn):
    redshift={}
    conn = psycopg2.connect("host=srv01050.soton.ac.uk dbname=frohmaier user=frohmaier password=rates")
    cur=conn.cursor()
    for sn in cclist:
        cur.execute("SELECT * from specinfo where ptfname=%s;", (str(sn[:-7]),))
        m=cur.fetchone()
       
        m=np.array(m)
        redshift[sn[:-7]] = m[5]
    cur.close()    
    return redshift[sn[:-7]]



#this function requires the model name, z and t0 values so that the values returned are specific to each data set
def modellimit(source_name, x, t0):
    source=sncosmo.get_source(source_name) 
    model=sncosmo.Model(source=source)
    
    #Creating a flux array of the model lightcurve
    timearray= np.arange(t0-120,t0+500, 0.01)
    model.set(z=spec_z(x),t0=t0)
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
        
        #adding zpsys and filter columns
        ab=np.zeros(len(hml[:,1]), dtype='|S2')
        for i in range(len(ab)):
            ab[i]='ab'
        hml=np.column_stack((hml,ab))   
        band=np.zeros(len(hml[:,1]), dtype='|S6')
        for i in range(len(band)):
            band[i]='ptf48r'
        hml=np.column_stack((hml,band))
        #fit to model
        z0 = float(spec_z(x[:-7]))
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(z0,z0 + 0.0001)}, nburn=10000, nsamples=50000)
        
        #The following excludes data points not in the range of the model and data sets with fewer than 4 data points
        limit = modellimit(source_name, x[:-7], res.parameters[1])
        hml2 = []
        for j in range(len(hml[:,1])):
            datapoint = hml [:,1][j]
            if (res.parameters[1]- limit[0])< float(datapoint) < (res.parameters[1]+limit[1]):
                hml2.append(hml[j]) 
        hml2 = np.array(hml2)
        if len(hml2)>3:
            
            return finalfitcurve(x, source_name, hml2) 
               
        
    except ValueError:
        print 'error'  
         
def finalfitcurve(x, source_name, hml):
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 
        
        #fit to model
        z0 = float(spec_z(x[:-7]))
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(z0, z0+0.0001)}, nburn=10000, nsamples=50000)
        #The following excludes data with not enough distribution..
        ##Ensuring at least 2 data points on either side of peak.
        limit = modellimit(source_name, x[:-7], res.parameters[1])
        j=0; k=0
        for i in hml[:,1]:
            if float(i)>float(limit[2]):
                j+=1
            else:
                k+=1
        if j>=2 and k>=2:          
            return np.array([float(res.chisq/res.ndof), hml[:,0][0], source_name]), hml #returns reduced chi squared, sn name and model name. 
           
    except ValueError:
        print 'error'

def showcurve(sn, source_name, hml):
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 

        #fit to model
        z0 = float(spec_z(sn[:-7]))
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(z0, z0+0.001)}, nburn=10000, nsamples=50000)

        #plot model
        sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name + '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=10)
        plt.show()

        #print 'Parameters:''z:',float(res.parameters[0]), float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])'
        print 'Done:', hml[:,0][0], 'z:',float(res.parameters[0]), 'Reduced chi^2:', res.chisq/res.ndof,#'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof

    except ValueError:
        print sn, source_name, 'cannot be plotted'

#creating list of available supernova data
datadir = '/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'

cclist = []
for i in os.walk(datadir):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)

print cclist
testlist = cclist[5:10]
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
        #find model with best fit determined by reduced chi squared
        minredchisq = min(redchisq)
        bestfit = fitted[redchisq.index(minredchisq)]
        
        #display plot
        
        #showcurve(sn, bestfit[0][2],bestfit[1])

