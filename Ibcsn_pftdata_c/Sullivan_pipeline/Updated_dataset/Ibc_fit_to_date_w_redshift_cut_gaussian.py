import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import psycopg2
from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time
from scipy import optimize
from astroML.plotting import hist

def averagedata(hml):
    duplicates = []
    totallen = len(hml[:,0])
    rownum = range(totallen)
    for i in rownum:
        lis = []
        for j in range(totallen):
            if i < j:
                if abs(float(hml[:,0][i])- float(hml[:,0][j])) < 0.25:
                    rownum.remove(j)
                    lis.append(j)
        lis.append(i)
        duplicates.append(lis)
    averagedata = []
    lenrow = len(hml[0])
    for lis in duplicates:
        alldata = [[] for z in range(lenrow)]
        for j in range(lenrow):
            alldata[j] = [hml[row][j] for row in lis]
            
        averagerow = [ z for z in range(lenrow)]
        counts = [1]
        errors = [2]
        others = [0, 3]
        for y in counts:
            datapoints = []
            for m in range(len(alldata[y])):
                datapoint = [alldata[y][m], alldata[y+1][m]]
                datapoints.append(datapoint) 
            denominator = np.sum([(1/np.square(float(x[1]))) for x in datapoints]) 
            numerator = np.sum([float(datapoint[0]) * (1/np.square(float(datapoint[1]))) for datapoint in datapoints])
            averagevalue = numerator/denominator
            averagerow[y]=averagevalue
            
        for y in errors:
            averagevalue = np.sqrt(np.sum([np.square((1/float(x[1])/denominator)) for x in datapoints ]))
            averagerow[y]=averagevalue
            
        for y in others:
            averagevalue = np.average([float(x) for x in alldata[y]])
            averagerow[y]=averagevalue
        averagedata.append(averagerow)
    hml = np.array(averagedata, dtype=str)
    hml2 = []
    for i in hml:
        if float(i[2])<1000:
            hml2.append(i)
        else:
            pass
    hml2 = np.array(hml2)
    return hml2
    
    
def modellimit(z0):
    source=sncosmo.get_source('nugent-sn1bc') 
    
    model=sncosmo.Model(source=source)
    t0 = 0
    #Creating a flux array of the model lightcurve
    timearray= np.arange(t0-120,t0+500, 0.01)
    model.set(z=z0,t0=t0)
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
            return timearray[list(fluxlc).index(max(fluxlc))]


def showcurve(sn, hml):

    
    source=sncosmo.get_source('nugent-sn1bc') 
    model=sncosmo.Model(source=source)
    z0 = float(redshifts[sn])
    model.set(z=z0)
    
    numofrows = len(hml[:,1])
    
    #adding flux, fluxerr, mag and magerror columns
    flux = [str(float(hml[i][1])*10**((-27-21.49)/2.5)) for i in range(numofrows)]
    hml=np.column_stack((hml,flux))
    fluxerr = [str(float(hml[i][2])*10**((-27-21.49)/2.5)) for i in range(numofrows)]
    hml=np.column_stack((hml,fluxerr))
    #adding zpsys and filter columns
    ab = ['ab' for i in range(numofrows)]
    hml=np.column_stack((hml,ab))  
    band = ['ptf48r' for i in range(numofrows)]
    hml=np.column_stack((hml,band))
    ''''
    'counts to flux eq',     hml_dat[8][1]*10**((-27-21.49)/2.5)
    'flux to mag eq',     -2.5 *np.log10(flux)-21.49
    'counts to mag eq',     -2.5 *np.log10(hml_dat[8][1]) +27
    '''
    #plt.errorbar([float(i) for i in hml[:,0]], [float(i) for i in hml[:,1]], yerr = [float(i) for i in hml[:,2]], color = 'blue', fmt = 'o')
    #plt.scatter(hml[:,0], hml[:,4], color = 'green')
    
    #plt.show()
    hml_dat=astropy.table.Table(data=hml, names=('mjd', 'counts', 'dcounts', 'zp', 'flux', 'fluxerr', 'zpsys', 'filter'), dtype=('float','float','float','float','float','float','str', 'str'))
    #print hml_dat
    
    mjdlimits = {'09': (55197-365-30, 55197),'10': (55197-30, 55562), '11': (55562-30, 55927), '12': (55927-30, 56293)}
    res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['t0','amplitude'], nburn=10000, nsamples=50000, bounds={'t0':mjdlimits[sn[:2]]})
    #res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['t0','amplitude'], nburn=10000, nsamples=50000, bounds={'t0':(55830,55960)})
    
   # sncosmo.plot_lc(hml_dat, model=fitted_model, model_label=source.name, errors=res.errors, color='blue', figtext=str(sn+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=18)
  #  plt.show()
    
    
    
    hml2 = []
    #print hml
    #print res.parameters[1]
    j = 0; k= 0
    start = float(res.parameters[1])
    peaktime = modellimit(z0)
    for row in hml:
        datestamp = row[0]
        if peaktime + (40.0)*(1+z0) > float(datestamp)-start > peaktime - (20.0)*(1+z0) :
            hml2.append(row)
        if peaktime + (40.0)*(1+z0) > float(datestamp)-start > peaktime:
            j+=1
        if peaktime > float(datestamp)-start> peaktime-(20.0)*(1+z0):
            k+=1  
            
    if j > 0 and k > 0:
        hml2= np.array(hml2)

        ## print 'after', len(hml2)
        hml_dat2=astropy.table.Table(data=hml2, names=('mjd', 'counts', 'dcounts', 'zp', 'flux', 'fluxerr', 'zpsys', 'filter'), dtype=('float','float','float','float','float','float','str', 'str'))
        
        res2, fitted_model2=sncosmo.fit_lc(hml_dat2, model, ['t0','amplitude'], nburn=10000, nsamples=50000, bounds={'t0':mjdlimits[sn[:2]]})

        #plot model
        redchisqaured = res2.chisq/res2.ndof

        #sncosmo.plot_lc(hml_dat2, model=fitted_model2, model_label=source.name, zp = 27.00, errors=res2.errors, color='blue', figtext=str(sn+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res2.chisq/res2.ndof), xfigsize=18)
       # plt.show()
        
        ##print 'Done:', hml2[:,0][0], 'z:',float(res2.parameters[0]), 'Reduced chi^2:', res2.chisq/res2.ndof #'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof
        return  redchisqaured, res2.parameters, z0, hml2 #, res.errors['amplitude']]
    else:
        print sn, 'didn\'t make the cut'
        #return res.parameters, np.mean(zp_new), res.chisq/res.ndof, z0, hml
# except ValueError:
    #pass
    
    #print sn, 'cannot be plotted'





completesnlist = ['09awk', '09dfk', '09dh', '09dha', '09dzt', '09fsr', '09iqd', '09ps', '09q', '09sk', '09ut', '10aavz', '10abck', '10acbu', '10acff', '10acgq', '10bhu', '10bip', '10bzf', '10ciw', '10cs', '10eqi', '10fbv', '10feq', '10fia', '10fmx', '10gvb', '10hfe', '10hie', '10inj', '10iue', '10kui', '10lbo', '10ood', '10osn', '10pbi', '10qif', '10qqd', '10qts', '10svt', '10tqi', '10tqv', '10vgv', '10vnv', '10wal', '10xem', '10xik', '10xjr', '10yow', '10ysd', '10zcn', '11awe', '11bli', '11bov', '11cmh', '11gcj', '11hyg', '11ilr', '11img', '11ixk', '11izq', '11jgj', '11kaa', '11klg', '11lbm', '11lmn', '11mnb', '11mwk', '11pnq', '11qcj', '11qiq', '11rfh', '11rka', '12as', '12bwq', '12cde', '12cjy', '12dcp', '12dtf', '12eaw', '12elh', '12fes', '12fgw', '12gdy', '12gps', '12grr', '12gty', '12gvr', '12gzk', '12hni', '12hvv', '12jxd', '12ktu', '12ldy', '12lpo', '12ltw', '12lvt', '12mfx']
types = {'10osn': 'SN Ic', '12gty': 'SN Ic', '10wal': 'SN Ic', '12mfx': 'SN Ib', '09ut': 'SN Ib/c', '10vgv': 'SN Ic-BL', '10acff': 'SN Ib', '09dh': 'SN Ic', '10xem': 'SN Ic-BL', '10aavz': 'SN Ic-BL', '11pnq': 'SN Ib/c', '12ktu': 'SN Ic', '12fgw': 'SN Ic', '10qts': 'SN Ic-BL', '11klg': 'SN Ic', '12gdy': 'SN Ic-BL', '09awk': 'SN Ib', '09q': 'SN Ic', '11cmh': 'SN Ic-BL', '12cde': 'SN Ib/c', '10hfe': 'SN Ic', '11ixk': 'SN Ic', '11bov': 'SN Ic', '10fia': 'SN Ib', '10cs': 'SN Ic-BL', '10qqd': 'SN Ic', '10bip': 'SN Ic', '10ysd': 'SN Ic-BL', '10kui': 'SN Ib', '12lpo': 'SN Ib/c', '12elh': 'SN Ic', '12gvr': 'SN Ib/c', '10fbv': 'SN Ib', '09dfk': 'SN Ib', '11img': 'SN Ic-BL', '11hyg': 'SN Ic', '10zcn': 'SN Ic', '12ldy': 'SN Ibn', '11izq': 'SN Ib', '11lmn': 'SN Ib/c', '09iqd': 'SN Ic', '11kaa': 'SN Ib', '09sk': 'SN Ic-BL', '10acgq': 'SN Ib', '11gcj': 'SN Ib/c', '12dcp': 'SN Ic', '10svt': 'SN Ic', '10vnv': 'SN Ib', '12eaw': 'SN Ib', '12as': 'SN Ic-BL', '12hni': 'SN Ic', '12gps': 'SN Ib', '10bhu': 'SN Ic', '09ps': 'SN Ic', '10abck': 'SN Ib', '11bli': 'SN Ib/c', '10pbi': 'SN Ib', '09fsr': 'SN Ib', '10hie': 'SN Ic', '10xik': 'SN Ic', '10inj': 'SN Ib', '10lbo': 'SN Ic', '10yow': 'SN Ic', '09dha': 'SN Ib', '12dtf': 'SN Ic', '12jxd': 'SN Ic', '10acbu': 'SN Ic', '10ciw': 'SN Ic-BL', '11lbm': 'SN Ic-BL', '11rka': 'SN Ic', '10feq': 'SN Ib', '10bzf': 'SN Ic-BL', '10gvb': 'SN Ic-BL', '11ilr': 'SN Ib', '12ltw': 'SN Ib', '10xjr': 'SN Ib', '10fmx': 'SN Ic', '10iue': 'SN Ic-BL', '10qif': 'SN Ib', '12fes': 'SN Ib', '12gzk': 'SN Ic', '12grr': 'SN Ic-BL', '11qcj': 'SN Ibn', '11rfh': 'SN Ibn', '11jgj': 'SN Ic', '10tqv': 'SN Ic', '12cjy': 'SN Ic', '09dzt': 'SN Ic', '11awe': 'SN Ib', '11mnb': 'SN Ic', '10eqi': 'SN Ib', '10ood': 'SN Ic', '10tqi': 'SN Ic', '12lvt': 'SN Ib', '12bwq': 'SN Ib', '11mwk': 'SN Ib/c', '11qiq': 'SN Ib', '12hvv': 'SN Ic'}
redshifts = {'10osn': 0.038, '12gty': 0.176, '10wal': 0.028803, '12mfx': 0.113, '09ut': 0.042, '10vgv': 0.015, '10acff': 0.06, '09dh': 0.07, '10xem': 0.0567, '10aavz': 0.062, '11pnq': 0.074, '12ktu': 0.031, '12fgw': 0.055, '10qts': 0.0907, '11klg': 0.026522, '12gdy': 0.157, '09awk': 0.062, '09q': 0.09, '11cmh': 0.1055, '12cde': 0.0125, '10hfe': 0.049, '11ixk': 0.021, '11bov': 0.022, '10fia': 0.039, '10cs': 0.166, '10qqd': 0.0807, '10bip': 0.051, '10ysd': 0.0963, '10kui': 0.021, '12lpo': 0.004486, '12elh': 0.116, '12gvr': 0.056, '10fbv': 0.056, '09dfk': 0.016, '11img': 0.158, '11hyg': 0.03, '10zcn': 0.02, '12ldy': 0.10598, '11izq': 0.062, '11lmn': 0.0904, '09iqd': 0.034, '11kaa': 0.04, '09sk': 0.0355, '10acgq': 0.1047, '11gcj': 0.148, '12dcp': 0.03093, '10svt': 0.031, '10vnv': 0.015, '12eaw': 0.028509, '12as': 0.033, '12hni': 0.107, '12gps': 0.016265, '10bhu': 0.036, '09ps': 0.1065, '10abck': 0.0143, '11bli': 0.03, '10pbi': 0.0475, '09fsr': 0.007942, '10hie': 0.067, '10xik': 0.071, '10inj': 0.0655, '10lbo': 0.053, '10yow': 0.02453, '09dha': 0.03, '12dtf': 0.061, '12jxd': 0.0254, '10acbu': 0.009877, '10ciw': 0.115, '11lbm': 0.039, '11rka': 0.0744, '10feq': 0.0279, '10bzf': 0.0498, '10gvb': 0.098, '11ilr': 0.0368, '12ltw': 0.06, '10xjr': 0.03, '10fmx': 0.047, '10iue': 0.1325, '10qif': 0.064, '12fes': 0.0359, '12gzk': 0.01377, '12grr': 0.16, '11qcj': 0.028, '11rfh': 0.0596, '11jgj': 0.04, '10tqv': 0.0795, '12cjy': 0.044, '09dzt': 0.0874, '11awe': 0.0553, '11mnb': 0.0603, '10eqi': 0.03, '10ood': 0.059, '10tqi': 0.0381, '12lvt': 0.011535, '12bwq': 0.04, '11mwk': 0.1215, '11qiq': 0.0317, '12hvv': 0.029}

datadir = '/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/Sullivan_pipeline/data_files/PTF'

tempcompletesnlist = ['09dfk', '09dh', '09dha','09dzt','09fsr','09q', '09sk', '09ut','10aavz', '10acbu', '10acff', '10acgq', '10bhu', '10bip', '10bzf', '10ciw', '10eqi', '10fbv', '10feq', '10fmx', '10gvb', '10hfe', '10hie', '10inj', '10iue', '10kui', '10lbo', '10osn', '10qif', '10qqd', '10qts', '10svt', '10tqi', '10tqv', '10vnv', '10xem', '10xik', '10xjr', '10ysd', '10zcn','11bli', '11bov', '11cmh', '11gcj', '11hyg', '11ilr', '11img', '11ixk', '11izq', '11jgj', '11kaa', '11klg', '11lbm', '11mnb', '11mwk', '11qcj', '11qiq', '11rka', '12as', '12cde', '12cjy', '12dcp', '12dtf', '12eaw', '12gdy', '12gzk', '12hni', '12jxd', '12ktu', '12lpo', '12mfx']
#tempcompletesnlistwithoutbadones
#snthatdontwork = #REMOVEDSN: '10fia','10wal','11rfh', '12bwq', #'12elh', '12fgw','12gty', '12lvt','11awe','11lmn', '11pnq', '12fes' ,'12grr','12gvr',  '12hvv','12ldy'
#completesnlist = tempcompletesnlistwithoutbadones
completesnlist = tempcompletesnlist
#completesnlist = snthatdontwork
#completesnlist = ['10fmx', '10xik', '11klg', '10vnv', '10xem', '10xjr', '10svt', '12eaw','11qiq']

source=sncosmo.get_source('nugent-sn1bc') 
model=sncosmo.Model(source=source)

ib = 0; ic = 0; jo = 1
redchib = []; redchic = []
redshiftvb = []; redshiftvc = []
totalred = []
peakabsmagvalueb = []; peakabsmagvaluec = []
for sn in completesnlist[:]:
    x = np.loadtxt(datadir + sn + '.extra_out_PTF48R', dtype= str,skiprows=16,  usecols = (0,1,2,3))# comments='#', delimiter='    ')
    dataset = []
    for i in x:
        if float(i[1]) > 0:
            dataset.append(i)
    hml = np.array(dataset)
    hml = averagedata(hml)
    
    try:
        y= showcurve(sn, hml)

        if y != None:
            redshiftlimit = 0.04
            redchilimit = 40
            if y[2] < redshiftlimit and y[0] < redchilimit:

                sntype = types[sn]
                if sntype[4] =='b':
                    ib+=1
                    redchib.append(y[0])
                else:
                    redchic.append(y[0])
                    ic+=1
                # redchi.append(y[2])
                #if y[0] > redchilimit:
                #   pass
                    #print sn, 'redchi: ', y[2]
                model.set(z=y[1][0],t0=y[1][1], amplitude=y[1][2]) 
                cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
                dispc = cosmo.luminosity_distance(y[2])*(10**(6))  #input redshift, will return in Mpc
                
                timearray = np.arange(y[1][1],y[1][1]+60,1)
                fluxlc = model.bandflux('ptf48r',timearray, zp = 27.0, zpsys='ab')
                
                
                obsmag =[(-21.49 - 2.5*np.log10(x)) for x in fluxlc]
                absmag =[(np.min(point) - 5*(np.log10(dispc.value) - 1)) for point in obsmag]
                
                fluxdatapoints = np.array(y[3][:,4], dtype = float)
                fluxdatapointsall = np.array([(float(i)*10**((-27-21.49)/2.5)) for i in hml[:,1]], dtype = float)
                
                fluxyerr = np.array(y[3][:,5], dtype = float)
                fluxyerrall = np.array([(float(i)*10**((-27-21.49)/2.5)) for i in hml[:,2]], dtype = float)
                
                
                obsdatapoints =np.array([-2.5 *np.log10(float(i)) +27 for i in y[3][:,1]], dtype = float)
                absdatapoints = [(np.min(point) - 5*(np.log10(dispc.value) - 1)) for point in obsdatapoints ]
                absmagvalue =np.min(absmag)
                
                obsmagyerr =[]
                for i in range(len(y[3][:,2])):
                    countserror = y[3][:,2]
                    counts = y[3][:,1]
                    obsmagyerr.append((2.5*float(countserror[i]))/(2.303*float(counts[i])))
                
                '''
                
                plt.subplot(3,3,jo)
                plt.plot(timearray, fluxlc, color= 'blue')
                plt.errorbar(np.array(hml[:,0], dtype = float), fluxdatapointsall , yerr=fluxyerrall, fmt = 'o', color= 'green')
                plt.errorbar(np.array(y[3][:,0], dtype = float), fluxdatapoints , yerr=fluxyerr, fmt = 'o', color= 'blue')
                plt.xlabel('Time /days')
                plt.ylabel('Flux') #Absolute Magnitude')
                plt.title(sn + ' absmag: ' + str(absmagvalue) + ' redchisq: ' + str(y[0]))
                plt.tight_layout()
                plt.gca().invert_xaxis()
                
                jo+=1
                '''
    
                print sn, 'abs peak mag', absmagvalue,
                
                if y[2] < redshiftlimit and y[0] < redchilimit:
                    totalred.append(y[2])
                    if sntype[4] =='b':
                        peakabsmagvalueb.append(absmagvalue)
                        redshiftvb.append(y[2])
                    else:
                        peakabsmagvaluec.append(absmagvalue)
                        redshiftvc.append(y[2])
    
                    
     
    except: 
         print sn, 'better data needed'           
                
            

#print jo
#plt.show()


#print j






#peakabsmagvalue = [peakabsmagvalueb,peakabsmagvaluec] 
#print peakabsmagvalueb
#hist(peakabsmagvalueb, bins = 'knuth', label = str(ib) + ' Ib datapoints', color = 'blue', histtype='stepfilled', alpha=0.2)#, stacked=True)


#hist(peakabsmagvaluec, bins = 'knuth', label = str(ic) + ' Ic datapoints', color ='green', histtype='stepfilled', alpha=0.2, des)#, stacked=True)



#plotting best fit gaussian
plt.subplot(221)  

result = hist(peakabsmagvalueb, bins = 'knuth', label = str(ib) + ' Ib datapoints', color = 'blue', histtype='stepfilled', alpha=0.2)#, stacked=True)

mean = np.mean(peakabsmagvalueb)
variance = np.var(peakabsmagvalueb)
sigma = np.sqrt(variance)
x = np.linspace(-21, -13,100)
dx = result[1][1] - result[1][0]
scale = len(peakabsmagvalueb)*dx

plt.plot(x, mlab.normpdf(x,mean,sigma)*scale, label = 'Best fit, mean: ' +str(round(mean, 3))+' sigma: ' + str(round(sigma, 3)))
 
print '\n','\n', 'Ib mean:', mean, 'sigma:', sigma

result = hist(peakabsmagvaluec, bins = 'knuth', label = str(ic) + ' Ic datapoints', color ='green', histtype='stepfilled', alpha=0.2)#, stacked=True)
#plt.xlim(-21, -13)
mean = np.mean(peakabsmagvaluec)
variance = np.var(peakabsmagvaluec)
sigma = np.sqrt(variance)
x = np.linspace(-21, -13,100)
dx = result[1][1] - result[1][0]
scale = len(peakabsmagvaluec)*dx



plt.plot(x, mlab.normpdf(x,mean,sigma)*scale, label = 'Best fit, mean: ' +str(round(mean, 3))+' sigma: ' + str(round(sigma, 3)))
print 'Ic mean:', mean, 'sigma:', sigma

plt.title("Peak Absolute Magnitude Histogram of TypeIb and TypeIc using the nugent-sn1bc model")
plt.xlabel("Absolute Magnitude")
plt.ylabel("Frequency")
plt.legend()



#Plotting Li gaussian
plt.subplot(223)  

result = hist(peakabsmagvalueb, bins = 'knuth', label = str(ib) + ' Ib datapoints', color = 'blue', histtype='stepfilled', alpha=0.2)#, stacked=True)
mean = -17.01
sigma = 0.41
x = np.linspace(-21, -13,100)
dx = result[1][1] - result[1][0]
scale = len(peakabsmagvalueb)*dx

plt.plot(x, mlab.normpdf(x,mean,sigma)*scale, label = 'mean: ' +str(round(mean, 3))+' sigma: ' + str(round(sigma, 3)))
 

result = hist(peakabsmagvaluec, bins = 'knuth', label = str(ic) + ' Ic datapoints', color ='green', histtype='stepfilled', alpha=0.2)#, stacked=True)
#plt.xlim(-21, -13)
mean = -16.04
sigma = 1.28
x = np.linspace(-21, -13,100)
dx = result[1][1] - result[1][0]
scale = len(peakabsmagvaluec)*dx

plt.plot(x, mlab.normpdf(x,mean,sigma)*scale, label =  'mean: ' +str(round(mean, 3))+' sigma: ' + str(round(sigma, 3)))



#hist(peakabsmagvaluec, bins = 6, label = str(ic) + ' Ic datapoints', color ='green' histtype='stepfilled', alpha=0.2)#, stacked=True)

#plt.hist(peakabsmagvaluec, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green',stacked=True)
plt.title("Peak Absolute Magnitude Histogram of TypeIb and TypeIc using the nugent-sn1bc model")
plt.xlabel("Absolute Magnitude")
plt.ylabel("Frequency")
plt.legend()




plt.subplot(222)            
redchi = [redchib,redchic]
hist(redchib, bins = 'knuth', label = str(ib) + ' Ib datapoints', color = 'blue', histtype='stepfilled', alpha=0.2)#, stacked=True)
hist(redchic, bins = 'knuth', label = str(ic) + ' Ic datapoints', color ='green' ,histtype='stepfilled', alpha=0.2)
#hist(redchi, bins = 10, label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'],histtype='stepfilled', alpha=0.2)#, stacked=True)
#plt.hist(redchic, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green')

plt.title("Reduced Chi Squared Distribution of TypeIb/c using the nugent-sn1bc model")
plt.xlabel("Reduced Chi Squared")
plt.ylabel("Frequency")
plt.legend()

plt.subplot(224)        


allredshifts = []
for i in redshifts:
    allredshifts.append( redshifts[i])
    
        
redshifts = [redshiftvb,redshiftvc]
totalredshifts = [allredshifts , ]
#hist(redshiftvb, bins = 'knuth', label = str(ib) + ' Ib datapoints' , color = 'blue', histtype='stepfilled', alpha=0.2)    #, stacked=True)
hist(allredshifts, bins = np.arange(0.00, 0.18,0.01) , label = str(len(allredshifts)) + ' sn available' , color = 'blue', histtype='stepfilled', alpha=0.2)    #, stacked=True)
hist(totalred, bins = np.arange(0.02, 0.1,0.01), label = str(ic +ib) + ' sn in sample used' , color = 'black', histtype='step')    #, stacked=True)

#hist(redshiftvc, bins = 'knuth', label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'])    #, stacked=True)

#plt.hist(redchic, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green')

plt.title("Redshift Distribution of TypeIb/c using the nugent-sn1bc model")
plt.xlabel("Redshift")
plt.ylabel("Frequency")
plt.legend()






plt.show()
