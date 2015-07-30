import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt

source=sncosmo.get_source('nugent-sn1bc',version='1.1') 
model=sncosmo.Model(source=source) 

hml=open('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt')  
fil = np.genfromtxt('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt', delimiter = '  ')#, dtype=(str, float, float,float, float,float, float,float,str,str,str,str))

lines = [line.split(' ') for line in hml if len(line)>2]
types = {}

for line in lines:
    for i in range(len(line)):
        if line[i]=='SN':
            types[line[0]] = line[i+1].strip(',')


def fitcurve(x):
    try:
        hml=np.load('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/' + x)    
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
    
        if types[hml[:,0][0]][1] == 'I':
            pass
        else:
            ##fitting model
            #print model
            res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['z','t0','amplitude'], bounds={'z':(0.005,0.35)}, nburn=10000, nsamples=50000)
            if len(hml[:,1])>3:

                #print res.errors
                sncosmo.plot_lc(hml_dat, model=fitted_model, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+' Type'+types[hml[:,0][0]]+'\n'+'Model name: '+ source.name),
                                                xfigsize=10)
                plt.show()
            #print #'### Parameters ###'
            #print 'SN',str(hml[:,0][0]), 'z:',float(res.parameters[0])# float(res.errors['z']), float(res.parameters[1]), float(res.errors['t0']),float(res.parameters[2]), float(res.errors['x0']),  float(res.parameters[3]), float(res.errors['x1']), float(res.parameters[4]), float(res.errors['c']), float(hml[:,8][0]), float(hml[:,9][0])
                print 'Done:', hml[:,0][0], 'z:',float(res.parameters[0]), 'Reduced chi^2:', res.chisq/res.ndof,'Data points:', len(hml[:,1]),'\n' #'Dof:', res.ndof

            else:
                pass
            #print 'Done:', hml[:,0][0], 'Not Ib/c probably'
            
    except ValueError:
        print 'Value Error' 
        


cclist = []
for i in os.walk('/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                fitcurve(k)
                
                cclist.append(k)
                
print len(cclist), 'supernovae in list'
#fitcurve('10hgi_LC.npy')