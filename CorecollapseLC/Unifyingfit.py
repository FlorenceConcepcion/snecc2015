import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt



x=open('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt')  
fil = np.genfromtxt('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt', delimiter = '  ')#, dtype=(str, float, float,float, float,float, float,float,str,str,str,str))

lines = [line.split(' ') for line in x if len(line)>2]
types = {}

for line in lines:
    for i in range(len(line)):
        if line[i]=='SN':
            types[line[0]] = line[i+1].strip(',')



def fitcurve(x, source_name):
    hml=np.load('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/' + x)    
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 
       # print len(hml[:,1]), 'data points'
       
        #adding zpsys and filter columns
        ab=np.zeros(len(hml[:,1]), dtype='|S2')
        for i in range(len(ab)):
            ab[i]='ab'
        hml=np.column_stack((hml,ab))   
        band=np.zeros(len(hml[:,1]), dtype='|S6')
        for i in range(len(band)):
            band[i]='ptf48r'
        hml=np.column_stack((hml,band))
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))

        ###fitting model
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(0.005,0.35)}, nburn=10000, nsamples=50000)
        #The following excludes data with fewer than 4 data points and not enough distribution.. 
        j=0; k=0
        for i in hml[:,1]:
            if float(i)>float(res.parameters[1]):
                j+=1
            else:
                k+=1
        # print hml[:,0][0],k,j
        if j>=2 and k>=2:            
            if len(hml[:,1])>3:                                
                #sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name), xfigsize=10)
                #plt.close()                
                
                return np.array([float(res.chisq/res.ndof), hml[:,0][0], source_name]) #returns reduced chi squared, sn name and model name. 
            else:
                pass
        else:
            pass
                
    except ValueError:
        pass
  #      print hml[:,0][0], 'Value Error with model: '+ source_name
        



def showcurve(sn, source_name):
    try:
        hml=np.load('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/' + sn)    
    
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source) 
    # print len(hml[:,1]), 'data points'
    
        #adding zpsys and filter columns
        ab=np.zeros(len(hml[:,1]), dtype='|S2')
        for i in range(len(ab)):
            ab[i]='ab'
        hml=np.column_stack((hml,ab))   
        band=np.zeros(len(hml[:,1]), dtype='|S6')
        for i in range(len(band)):
            band[i]='ptf48r'
        hml=np.column_stack((hml,band))
        hml_dat=astropy.table.Table(data=hml, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
        
        ###fitting model
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(0.005,0.35)}, nburn=10000, nsamples=50000)
            #print 'peak', float(res.parameters[1])
        print res                              
        sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name), xfigsize=10)
        #plt.clfig()
        plt.show()

        #print #'### Parameters ###'
        #print 'SN',str(hml[:,0][0]), 'z:',float(res.parameters[0])# float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])
        print 'Done:', hml[:,0][0], 'z:',float(res.parameters[0]), 'Reduced chi^2:', res.chisq/res.ndof,#'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof

    except ValueError:
        sn, source_name, 'cannot be plotted'

cclist = []
for i in os.walk('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
                
sourcelist= ['nugent-sn1a', 'nugent-sn91t', 'nugent-sn91bg', 'nugent-sn1bc', 'nugent-hyper', 'nugent-sn2p', 'nugent-sn2l', 'nugent-sn2n', 's11-2004hx', 's11-2005lc', 's11-2005hl', 's11-2005hm', 's11-2005gi', 's11-2006fo', 's11-2006jo', 's11-2006jl', 'hsiao', 'hsiao', 'hsiao', 'hsiao-subsampled', 'salt2', 'salt2', 'salt2', 'salt2', 'salt2-extended', 'snf-2011fe', 'snana-2004fe', 'snana-2004gq', 'snana-sdss004012', 'snana-2006fo', 'snana-sdss014475', 'snana-2006lc', 'snana-2007ms', 'snana-04d1la', 'snana-04d4jv', 'snana-2004gv', 'snana-2006ep', 'snana-2007y', 'snana-2004ib', 'snana-2005hm', 'snana-2007nc', 'snana-2004hx', 'snana-2005gi', 'snana-2006gq', 'snana-2006kn', 'snana-2006jl', 'snana-2006iw', 'snana-2006kv', 'snana-2006ns', 'snana-2007iz', 'snana-2007nr', 'snana-2007kw', 'snana-2007ky', 'snana-2007lj', 'snana-2007lb', 'snana-2007ll', 'snana-2007nw', 'snana-2007ld', 'snana-2007md', 'snana-2007lz', 'snana-2007lx', 'snana-2007og', 'snana-2007ny', 'snana-2007nv', 'snana-2007pg', 'snana-2006ez', 'snana-2006ix', 'whalen-z15b', 'whalen-z15d', 'whalen-z15g', 'whalen-z25b', 'whalen-z25d', 'whalen-z25g'] 
#, 'snana-2006jo' 'whalen-z40b', 'whalen-z40g' removed

testlist = cclist[15:20]
testlist = ['10svt_LC.npy']
for sn in testlist:
    fitted=[]; redchisq=[]
    for source in sourcelist:
        if fitcurve(sn, source)!= None:
            fitted.append(fitcurve(sn, source))
            redchisq.append(float(fitcurve(sn, source)[0]))        
        else:
            pass
    if len(fitted)!=0:
        minredchisq = min(redchisq)
        bestfit = fitted[redchisq.index(minredchisq)]
        print sn[:-4], bestfit[2]
        showcurve(sn, bestfit[2])
    else:
        pass

