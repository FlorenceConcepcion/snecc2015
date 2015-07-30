import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
import psycopg2
from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time



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
    return hml

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
    for i in range(1):
        
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
        #plt.scatter(hml[:,0], hml[:,1], color = 'blue')
        #plt.scatter(hml[:,0], hml[:,4], color = 'green')
        
        #plt.show()
        hml_dat=astropy.table.Table(data=hml, names=('mjd', 'counts', 'dcounts', 'zp', 'flux', 'fluxerr', 'zpsys', 'filter'), dtype=('float','float','float','float','float','float','str', 'str'))
        print hml_dat
        
        res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['t0','amplitude'], nburn=10000, nsamples=50000)
        ##sncosmo.plot_lc(hml_dat, model=fitted_model, model_label=source.name, errors=res.errors, color='blue', figtext=str(sn+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=18)
        #plt.show()
        hml2 = []
        #print hml
        #print res.parameters[1]
        j = 0; k= 0
        start = float(res.parameters[1])
        peaktime = modellimit(z0)
        for row in hml:
            datestamp = row[0]
            if peaktime + (40.0)*(1+z0) > float(datestamp)-start > peaktime - (15.0)*(1+z0) :
                hml2.append(row)
            if peaktime + (40.0)*(1+z0) > float(datestamp)-start > peaktime:
                j+=1
            if peaktime > float(datestamp)-start> peaktime-(15.0)*(1+z0):
                k+=1  
        if j > 1 and k>1:
            hml2= np.array(hml2)

            ##print 'after', len(hml2)
            hml_dat2=astropy.table.Table(data=hml2, names=('mjd', 'counts', 'dcounts', 'zp', 'flux', 'fluxerr', 'zpsys', 'filter'), dtype=('float','float','float','float','float','float','str', 'str'))
            res2, fitted_model2=sncosmo.fit_lc(hml_dat2, model, ['t0','amplitude'], nburn=10000, nsamples=50000)
            
            #plot model
            redchisqaured = res2.chisq/res2.ndof

            #sncosmo.plot_lc(hml_dat2, model=fitted_model2, model_label=source.name, zp = 27.00, errors=res2.errors, color='blue', figtext=str(sn+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res2.chisq/res2.ndof), xfigsize=18)
            #plt.show()
            
            ##print 'Done:', hml2[:,0][0], 'z:',float(res2.parameters[0]), 'Reduced chi^2:', res2.chisq/res2.ndof #'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof
            return  redchisqaured, res2.parameters, z0, hml2 #, res.errors['amplitude']]
        else:
            print sn, 'need better data'
            #return res.parameters, np.mean(zp_new), res.chisq/res.ndof, z0, hml
   # except ValueError:
        #pass
        print sn, 'cannot be plotted'





completesnlist = ['10osn', '10ciw', '12gty', '10wal', '12mfx', '10vgv', '10acff', '10xem', '10bhu', '11pnq', '12ktu', '12fgw', '10qts', '11klg', '12gdy', '11cmh', '12cde', '10hfe', '11ixk', '11bov', '10vwg', '10ysd', '10cs', '10bip', '10fia', '11oni', '10acgq', '12lpo', '12elh', '12cjy', '10hcw', '12gvr', '10fbv', '11img', '10zcn', '12ldy', '11izq', '11lmn', '11kaa', '10kui', '12fsu', '11gcj', '12dcp', '10svt', '10vnv', '12eaw', '12as', '12hni', '12gps', '10aavz', '10qqd', '11bli', '10pbi', '11hyg', '10hie', '10xik', '10inj', '10lbo', '10yow', '10qif', '12jxd', '10acbu', '10uhf', '11lbm', '10uhn', '11rka', '10abck', '10feq', '10bzf', '10gvb', '11ilr', '12ltw', '10xjr', '10fmx', '10iue', '10hgi', '12dtf', '12fes', '12gzk', '12grr', '10xfh', '11qcj', '11rfh', '11jgj', '10bfz', '10tqv', '11all', '11awe', '11mnb', '10eqi', '10ood', '10tqi', '12lvt', '12bwq', '11dhf', '11mwk', '11qiq', '12hvv']
types = {'10osn': 'Ic', '10ciw': 'Ic-BL', '12gty': 'Ic', '10wal': 'Ic', '12mfx': 'Ib', '10vgv': 'Ic-BL', '10acff': 'Ib', '10xem': 'Ic-BL', '10aavz': 'Ic-BL', '11pnq': 'Ib/c', '12ktu': 'Ic', '10qqd': 'Ic', '10qts': 'Ic-BL', '11klg': 'Ic', '12gdy': 'Ic-BL', '11cmh': 'Ic-BL', '12cde': 'Ib/c', '10hfe': 'Ic', '11ixk': 'Ic', '11bov': 'Ic', '10vwg': 'Ic', '10ysd': 'Ic-BL', '10cs': 'Ic-BL', '10bip': 'Ic', '10fia': 'Ib', '11oni': 'Ib/c', '10kui': 'Ib', '12lpo': 'Ib/c', '10pbi': 'Ib', '11all': 'Ic', '10hcw': 'Ib', '12gvr': 'Ib/c', '10fbv': 'Ib', '11img': 'Ic-BL', '10zcn': 'Ic', '12ldy': 'Ibn', '11izq': 'Ib', '11lmn': 'Ib/c', '11kaa': 'Ib', '10acgq': 'Ib', '12fsu': 'Ic', '11gcj': 'Ib/c', '12dcp': 'Ibn', '10svt': 'Ic', '10vnv': 'Ib', '12eaw': 'Ib', '12as': 'Ic-BL', '12hni': 'Ic', '12gps': 'Ib', '10bhu': 'Ic', '12fgw': 'Ic', '11bli': 'Ib/c', '12elh': 'Ic', '11hyg': 'Ic', '10hie': 'Ic', '10xik': 'Ic', '10inj': 'Ib', '10lbo': 'Ic', '10yow': 'Ic', '10qif': 'Ib', '12jxd': 'Ic', '10acbu': 'Ic', '10uhf': 'Ic', '11lbm': 'Ic-BL', '10uhn': 'Ic', '11rka': 'Ic', '10abck': 'Ib', '10feq': 'Ib', '10bzf': 'Ic-BL', '10gvb': 'Ic-BL', '11ilr': 'Ib', '12ltw': 'Ib', '10xjr': 'Ib', '10fmx': 'Ic', '10iue': 'Ic-BL', '10hgi': 'Ib/c', '12dtf': 'Ic', '12fes': 'Ib', '12gzk': 'Ic', '12grr': 'Ic-BL', '10xfh': 'Ic', '11qcj': 'Ic', '11rfh': 'Ic', '11jgj': 'Ic', '10bfz': 'Ic-BL', '10tqv': 'Ic', '12cjy': 'Ic', '11awe': 'Ib', '11mnb': 'Ic', '10eqi': 'Ib', '10ood': 'Ic', '10tqi': 'Ic', '12lvt': 'Ib', '12bwq': 'Ib', '11dhf': 'Ib', '11mwk': 'Ib/c', '11qiq': 'Ib', '12hvv': 'Ic'}
datadir = '/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/Sullivan_pipeline/data_files/PTF'
redshifts = {'10osn': '0.038', '10ciw': '0.1311', '12gty': '0.176', '10wal': '0.028803', '12mfx': '0.113', '10vgv': '0.015', '10acff': '0.06', '10xem': '0.0567', '10aavz': '0.062', '11pnq': '0.074', '12ktu': '0.031', '10qqd': '0.0807', '10qts': '0.0907', '11klg': '0.026522', '12gdy': '0.157', '11cmh': '0.1', '12cde': '0.0125', '10hfe': '0.049', '11ixk': '0.021', '11bov': '0.022', '10vwg': '0.19', '10ysd': '0.0963', '10cs': '0.166', '10bip': '0.051', '10fia': '0.039', '11oni': '99999.0', '10kui': '0.021', '12lpo': '0.004486', '12elh': '0.116', '11all': '0.19', '10hcw': '0.011128', '12gvr': '0.056', '10fbv': '0.056', '11img': '0.158', '10zcn': '0.02', '12ldy': '0.10598', '11izq': '0.062', '11lmn': '0.0904', '11kaa': '0.04', '10acgq': '0.1047', '12fsu': '0.1755', '11gcj': '0.148', '12dcp': '0.03093', '10svt': '0.031', '10vnv': '0.015', '12eaw': '0.028509', '12as': '0.033', '12hni': '0.107', '12gps': '0.016265', '10bhu': '0.036', '12fgw': '0.055', '11bli': '0.03', '10pbi': '0.0475', '11hyg': '0.03', '10hie': '0.067', '10xik': '0.071', '10inj': '0.0655', '10lbo': '0.053', '10yow': '0.02453', '10qif': '0.064', '12jxd': '0.0254', '10acbu': '0.009877', '10uhf': '0.289', '11lbm': '0.039', '10uhn': '0.101', '11rka': '0.0744', '10abck': '0.0143', '10feq': '0.0279', '10bzf': '0.0498', '10gvb': '0.1', '11ilr': '0.0368', '12ltw': '0.06', '10xjr': '0.03', '10fmx': '0.047', '10iue': '0.1325', '10hgi': '0.0985', '12dtf': '0.061', '12fes': '0.0359', '12gzk': '0.01377', '12grr': '0.16', '10xfh': '0.0174', '11qcj': '0.028', '11rfh': '0.0596', '11jgj': '0.04', '10bfz': '0.169', '10tqv': '0.0795', '12cjy': '0.044', '11awe': '0.0553', '11mnb': '0.0603', '10eqi': '0.03', '10ood': '0.059', '10tqi': '0.0381', '12lvt': '0.011535', '12bwq': '0.04', '11dhf': '0.028', '11mwk': '0.1215', '11qiq': '0.0317', '12hvv': '0.029'}


tempcompletesnlist = ['10osn', '10ciw', '12gty', '10wal', '12mfx', '10xem', '10bhu', '11pnq', '12ktu', '12fgw', '10qts', '11klg', '11cmh', '12cde', '10hfe', '11ixk', '11bov', '10ysd', '10bip', '10fia', '11oni', '12lpo', '12elh', '12cjy', '10hcw', '12gvr', '10fbv', '11img', '10zcn', '11izq', '10kui', '12fsu', '11gcj', '12dcp', '10svt', '10vnv', '12eaw', '12hni', '10qqd', '11bli', '11hyg', '10hie', '10xik', '10inj', '10lbo', '10qif', '12jxd', '10uhf', '10uhn', '11rka', '10feq', '10bzf', '10gvb', '11ilr', '10xjr', '10fmx', '10iue', '10hgi', '12dtf', '12fes', '12gzk', '12grr', '11qcj', '11rfh', '11jgj', '10bfz', '10tqv', '11all', '11mnb', '10eqi', '10tqi', '12lvt', '12bwq', '11mwk', '11qiq', '12hvv']

completesnlist = tempcompletesnlist
#completesnlist = ['11pnq']
source=sncosmo.get_source('nugent-sn1bc') 
model=sncosmo.Model(source=source)

ib = 0; ic = 0; jo = 1
redchib = []; redchic = []
peakabsmagvalueb = []; peakabsmagvaluec = []
for sn in completesnlist[:10]:
    print sn
    x = np.loadtxt(datadir + sn + '.extra_out_PTF48R', dtype= str,skiprows=16,  usecols = (0,1,2,3))# comments='#', delimiter='    ')
    dataset = []
    for i in x:
        if float(i[1]) > 0:
            dataset.append(i)
    hml = np.array(dataset)
    hml = averagedata(hml)
    
    try:
        y= showcurve(sn, hml)
        print y

        if y != None:
            redchilimit = 2000
            if y[0] < redchilimit:
                sntype = types[sn]
                if sntype[1] =='b':
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
                
                
                        
                plt.subplot(3,4,jo)
                plt.plot(timearray, fluxlc, color= 'blue')
                plt.errorbar(np.array(hml[:,0], dtype = float), fluxdatapointsall , yerr=fluxyerrall, fmt = 'o', color= 'green')
                plt.errorbar(np.array(y[3][:,0], dtype = float), fluxdatapoints , yerr=fluxyerr, fmt = 'o', color= 'blue')
                plt.xlabel('Time /days')
                plt.ylabel('Flux') #Absolute Magnitude')
                plt.title(sn + ' absmag: ' + str(absmagvalue) + ' redchisq: ' + str(y[0]))
            # plt.tight_layout()
            # plt.gca().invert_yaxis()
                
                jo+=1
                
    
    
                if y[0] < redchilimit:
                    if sntype[1] =='b':
                        peakabsmagvalueb.append(absmagvalue)
                    else:
                        peakabsmagvaluec.append(absmagvalue)
    
                    print sn, 'abs peak mag', absmagvalue
     
    except: 
         pass           
                
            

#print jo
plt.show()



#print j

'''    
plt.subplot(211)  
peakabsmagvalue = [peakabsmagvalueb,peakabsmagvaluec] 
plt.hist(peakabsmagvalue, bins = 10, label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'], stacked=True)
#plt.hist(peakabsmagvaluec, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green',stacked=True)
plt.title("Peak Absolute Magnitude Histogram of TypeIb and TypeIc using the nugent-sn1bc model")
plt.xlabel("Absolute Magnitude")
plt.ylabel("Frequency")
plt.legend()

plt.subplot(212)            
redchi = [redchib,redchic]
plt.hist(redchi, bins = 10, label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'], stacked=True)
#plt.hist(redchic, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green')

plt.title("Reduced Chi Squared Distribution of TypeIb/c using the nugent-sn1bc model")
plt.xlabel("Reduced Chi Squared")
plt.ylabel("Frequency")
plt.legend()
plt.show()
'''

