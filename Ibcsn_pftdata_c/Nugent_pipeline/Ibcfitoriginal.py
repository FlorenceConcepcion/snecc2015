import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
import psycopg2
from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time

#creating a dictionary with supernova name and type of supernova  
'''
x=open('/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/SNlist.txt')  

lines = [line.split(' ') for line in x if len(line)>2]
types = {}

for line in lines:
    for i in range(len(line)):
        if line[i]=='SN':
            types[line[0][:5]] = line[i+1].strip(',')
'''
types = {'10osn': 'Ic', '10ciw': 'Ic-BL', '12gty': 'Ic', '10wal': 'Ic', '12mfx': 'Ib', '10vgv': 'Ic-BL', '10xem': 'Ic-BL', '10bhu': 'Ic', '11pnq': 'Ib/c', '10abc': 'Ib', '12ktu': 'Ic', '10qqd': 'Ic', '10qts': 'Ic-BL', '11klg': 'Ic', '12gdy': 'Ic-BL', '11cmh': 'Ic-BL', '12cde': 'Ib/c', '10hfe': 'Ic', '11ixk': 'Ic', '10aav': 'Ic-BL', '11bov': 'Ic', '10vwg': 'Ic', '10ysd': 'Ic-BL', '10cs': 'Ic-BL', '10bip': 'Ic', '10fia': 'Ib', '11oni': 'Ib/c', '10kui': 'Ib', '12lpo': 'Ib/c', '10pbi': 'Ib', '11all': 'Ic', '10hcw': 'Ib', '12gvr': 'Ib/c', '10fbv': 'Ib', '11img': 'Ic-BL', '10zcn': 'Ic', '12ldy': 'Ibn', '11izq': 'Ib', '11lmn': 'Ib/c', '11kaa': 'Ib', '12fsu': 'Ic', '11gcj': 'Ib/c', '10acf': 'Ib', '10acg': 'Ib', '10acb': 'Ic', '12dcp': 'Ibn', '10svt': 'Ic', '10vnv': 'Ib', '12eaw': 'Ib', '12as': 'Ic-BL', '12hni': 'Ic', '12gps': 'Ib', '12fgw': 'Ic', '11bli': 'Ib/c', '12elh': 'Ic', '11hyg': 'Ic', '10hie': 'Ic', '10xik': 'Ic', '10inj': 'Ib', '10lbo': 'Ic', '10yow': 'Ic', '10qif': 'Ib', '12jxd': 'Ic', '10uhf': 'Ic', '11lbm': 'Ic-BL', '10uhn': 'Ic', '11rka': 'Ic', '10feq': 'Ib', '10bzf': 'Ic-BL', '10gvb': 'Ic-BL', '11ilr': 'Ib', '12ltw': 'Ib', '10xjr': 'Ib', '10fmx': 'Ic', '10iue': 'Ic-BL', '10hgi': 'Ib/c', '12dtf': 'Ic', '12fes': 'Ib', '12gzk': 'Ic', '12grr': 'Ic-BL', '10xfh': 'Ic', '11qcj': 'Ic', '11rfh': 'Ic', '11jgj': 'Ic', '10bfz': 'Ic-BL', '10tqv': 'Ic', '12cjy': 'Ic', '11awe': 'Ib', '11mnb': 'Ic', '10eqi': 'Ib', '10ood': 'Ic', '10tqi': 'Ic', '12lvt': 'Ib', '12bwq': 'Ib', '11dhf': 'Ib', '11mwk': 'Ib/c', '11qiq': 'Ib', '12hvv': 'Ic'}

'''
#get spectroscopic redshift from database
def spec_z(sn):
    conn = psycopg2.connect("host=srv01050.soton.ac.uk dbname=frohmaier user=frohmaier password=rates")
    cur=conn.cursor()
    cur.execute("SELECT * from specinfo where ptfname=%s;", (str(sn[:-7]),))
    m=cur.fetchone()
    m=np.array(m)
    cur.close() 
    try:
        return m[5]
    except IndexError:
        return 'no data'
'''     
specs = {'10xik_LC.npy': '0.071', '11rka_LC.npy': '0.0744', '10iue_LC.npy': '0.1325', '11bov_LC.npy': '0.022', '10svt_LC.npy': '0.031', '12ktu_LC.npy': '0.031', '11qiq_LC.npy': '0.0317', '12gzk_LC.npy': '0.01377', '10xem_LC.npy': '0.0567', '10vgv_LC.npy': '0.015', '12fes_LC.npy': '0.0359', '12dcp_LC.npy': '0.03093', '12fsu_LC.npy': '0.1755', '12cjy_LC.npy': '0.044', '10hgi_LC.npy': '0.0985', '10lbo_LC.npy': '0.053', '11kaa_LC.npy': '0.04', '12elh_LC.npy': '0.116', '10inj_LC.npy': '0.0655', '12gdy_LC.npy': '0.157', '10feq_LC.npy': '0.0279', '12fgw_LC.npy': '0.055', '12eaw_LC.npy': '0.028509', '11img_LC.npy': '0.158', '10qts_LC.npy': '0.0907', '10tqi_LC.npy': '0.0381', '11ilr_LC.npy': '0.0368', '10qqd_LC.npy': '0.0807', '11izq_LC.npy': '0.062', '12cde_LC.npy': '0.0125', '10vnv_LC.npy': '0.015', '11cmh_LC.npy': '0.1', '12dtf_LC.npy': '0.061', '12hni_LC.npy': '0.107', '12hvv_LC.npy': '0.029', '10ciw_LC.npy': '0.1311', '10xfh_LC.npy': '0.0174', '11mnb_LC.npy': '0.0603'}
        
              
def modellimit(z0):
    source=sncosmo.get_source('nugent-sn1bc') 
    t0 = 0
    model=sncosmo.Model(source=source)
    
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
            
            

def showcurve(sn, hml, source_name):
    try:
        source=sncosmo.get_source(source_name) 
        model=sncosmo.Model(source=source)
        z0 = specs[sn]
        if z0 == 'no data':
            #pass
            print sn[:-7], 'cannot be plotted'
        else:
            z0 = float(z0) 
            model.set(z=z0)
            #print z0
        
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
            res, fitted_model=sncosmo.fit_lc(hml_dat, model, ['t0','amplitude'], nburn=10000, nsamples=50000)
            #sncosmo.plot_lc(hml_dat, model=fitted_model, model_label=source.name, errors=res.errors, color='blue', figtext=str(hml[:,0][0]+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res.chisq/res.ndof), xfigsize=18)
            #plt.show()
            zp_new = [hml_dat[i][7] for i in range(len(hml))]
            
            
            hml2 = []
            #print hml
            #print res.parameters[1]
            j = 0; k= 0
            for row in hml:
                datestamp = row[1]
                start = float(res.parameters[1])
                peaktime = modellimit(z0)  
                if (40.0)*(1+z0) > float(datestamp) - start > peaktime-10 :
                    hml2.append(row)
                if (40.0)*(1+z0) > float(datestamp) - start > peaktime:
                    j+=1
                if peaktime > float(datestamp) - start > peaktime-10:
                    k+=1                         
            if j > 1 and k>1:
                hml2= np.array(hml2)

                ##print 'after', len(hml2)
                hml_dat2=astropy.table.Table(data=hml2, names=('ptfname', 'time', 'magnitude', 'mag_err', 'flux', 'flux_err', 'zp_new', 'zp', 'ra', 'dec', 'zpsys', 'filter'), dtype=('str','float','float','float','float','float','float','float','float','float','str', 'str'))
                res2, fitted_model2=sncosmo.fit_lc(hml_dat2, model, ['t0','amplitude'], nburn=10000, nsamples=50000)
                
                #plot model
                redchisqaured = res2.chisq/res2.ndof
                if redchisqaured < 1:
                    pass
                #sncosmo.plot_lc(hml_dat2, model=fitted_model2, model_label=source.name, zp = np.mean(zp_new), errors=res2.errors, color='blue', figtext=str(hml[:,0][0]+'\n'+'Model name: nugent-sn1bc'+ '\n'+'Reduced Chi Squared: ')+ str(res2.chisq/res2.ndof), xfigsize=18)
                #plt.show()
                
                ##print 'Done:', hml2[:,0][0], 'z:',float(res2.parameters[0]), 'Reduced chi^2:', res2.chisq/res2.ndof #'Data points:', len(hml[:,1]),' Type'+types[hml[:,0][0]]+'Model name: '+ source.name,'\n' #'Dof:', res.ndof
                return res2.parameters, np.mean(zp_new), redchisqaured, z0, hml2 #, res.errors['amplitude']]
            else:
                pass
                #return res.parameters, np.mean(zp_new), res.chisq/res.ndof, z0, hml
    except ValueError:
        #pass
        print sn, 'cannot be plotted'
    

#creating list of available supernova data
datadir = '/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/'
'''
cclist = []
for i in os.walk(datadir):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)

testlist = cclist[:]
'''
testlist = ['12ltw_LC.npy', '12cjy_LC.npy', '10bzf_LC.npy', '11qiq_LC.npy', '12cde_LC.npy', '11mwk_LC.npy', '10hie_LC.npy', '10hcw_LC.npy', '12lpo_LC.npy', '10qts_LC.npy', '10zcn_LC.npy', '10aav_LC.npy', '11hyg_LC.npy', '10osn_LC.npy', '10fia_LC.npy', '11all_LC.npy', '12grr_LC.npy', '10vgv_LC.npy', '12eaw_LC.npy', '10hgi_LC.npy', '11bov_LC.npy', '11oni_LC.npy', '10ysd_LC.npy', '10ciw_LC.npy', '10bip_LC.npy', '10qqd_LC.npy', '10gvb_LC.npy', '10xik_LC.npy', '11klg_LC.npy', '11kaa_LC.npy', '12mfx_LC.npy', '10svt_LC.npy', '10uhn_LC.npy', '11jgj_LC.npy', '11bli_LC.npy', '11ixk_LC.npy', '12gvr_LC.npy', '12dcp_LC.npy', '10bhu_LC.npy', '10acf_LC.npy', '10vwg_LC.npy', '10wal_LC.npy', '11img_LC.npy', '12lvt_LC.npy', '12hvv_LC.npy', '10inj_LC.npy', '10kui_LC.npy', '10feq_LC.npy', '12fsu_LC.npy', '10tqi_LC.npy', '12bwq_LC.npy', '10hfe_LC.npy', '12jxd_LC.npy', '10yow_LC.npy', '11gcj_LC.npy', '10xjr_LC.npy', '12fgw_LC.npy', '10acg_LC.npy', '12fes_LC.npy', '10eqi_LC.npy', '11pnq_LC.npy', '11qcj_LC.npy', '11mnb_LC.npy', '10pbi_LC.npy', '11rfh_LC.npy', '11izq_LC.npy', '10tqv_LC.npy', '10ood_LC.npy', '12gps_LC.npy', '11ilr_LC.npy', '12gty_LC.npy', '10bfz_LC.npy', '10qif_LC.npy', '10vnv_LC.npy', '10xem_LC.npy', '10xfh_LC.npy', '12elh_LC.npy', '10lbo_LC.npy', '12hni_LC.npy', '11cmh_LC.npy', '10fbv_LC.npy', '12gdy_LC.npy', '10iue_LC.npy', '12ktu_LC.npy', '10abc_LC.npy', '10uhf_LC.npy', '12gzk_LC.npy', '10acb_LC.npy', '10fmx_LC.npy', '11dhf_LC.npy', '12dtf_LC.npy', '11rka_LC.npy']


testlist = ['10vnv_LC.npy','11qiq_LC.npy', '10qts_LC.npy']
ib=ic=0
reducedlist = ['12cjy_LC.npy', '11qiq_LC.npy', '12cde_LC.npy', '10qts_LC.npy', '10vgv_LC.npy', '12eaw_LC.npy', '10hgi_LC.npy', '11bov_LC.npy', '10ciw_LC.npy', '10qqd_LC.npy', '10xik_LC.npy', '11kaa_LC.npy', '10svt_LC.npy', '12dcp_LC.npy', '11img_LC.npy', '12hvv_LC.npy', '10inj_LC.npy', '10feq_LC.npy', '12fsu_LC.npy', '10tqi_LC.npy', '12fgw_LC.npy', '12fes_LC.npy', '11mnb_LC.npy', '11izq_LC.npy', '11ilr_LC.npy', '10vnv_LC.npy', '10xem_LC.npy', '10xfh_LC.npy', '12elh_LC.npy', '10lbo_LC.npy', '12hni_LC.npy', '11cmh_LC.npy', '12gdy_LC.npy', '10iue_LC.npy', '12ktu_LC.npy', '12gzk_LC.npy', '12dtf_LC.npy', '11rka_LC.npy']
testlist = reducedlist[:]
print len(testlist)
#datalistmjd = []
redchic = []; redchib = []
peakabsmagvalueb=[]; peakabsmagvaluec=[]
source_name = 'nugent-sn1bc'
source=sncosmo.get_source(source_name) 
model=sncosmo.Model(source=source)
jo=1
for sn in testlist:

    #averaging rows of data taken within the same night
    hml = np.load(datadir + sn)
    duplicates = []; dup = {}
    totallen = len(hml[:,1])
    rownum = range(len(hml[:,1])-1)

    for i in rownum:
        lis = []
        for j in range(len(hml[:,1])-1):
            if i < j:
                if abs(float(hml[:,1][i])- float(hml[:,1][j])) < 0.25:
                    rownum.remove(j)
                    lis.append(j)
                    #print abs(float(hml[:,1][i])- float(hml[:,1][j]))
        lis.append(i)
        duplicates.append(lis)
    averagedata = []
    for lis in duplicates:
        alldata = [[] for z in range(10)]
        for j in range(10):
            alldata[j] = [hml[row][j] for row in lis]
        averagerow = [hml[0][0] ]+ [ n for n in range(9)]
        fluxmag = [2, 4]
        errors = [3, 5]
        others = [1, 6,7,8,9]
        for y in fluxmag:
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
    
    hml = np.array(averagedata[:])
    #print hml
    
    #z= np.min([float(x) for x in hml[:,1]])
    #datalistmjd.append([sn[:-7], z-50])
    
    if len(hml) >4:
        y = showcurve(sn, hml, source_name)
        if y != None:
            redchilimit = 9
            if y[2] < redchilimit:
                
                sntype = types[sn[:-7]]
                if sntype[1] =='b':
                    #print sntype
                    ib+=1
                    redchib.append(y[2])
                else:
                    redchic.append(y[2])
                    ic+=1
               # redchi.append(y[2])
            if y[2] >redchilimit:
                pass
                #print sn, 'redchi: ', y[2]
            model.set(z=y[0][0],t0=y[0][1], amplitude=y[0][2]) 
            cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
            dispc = cosmo.luminosity_distance(y[3])*(10**(6))  #input redshift, will return in Mpc
            
            timearray = np.arange(y[0][1],y[0][1]+60,1)
            fluxlc = model.bandflux('ptf48r',timearray, zp = y[1],zpsys='ab' )
            obsmag =[(float(y[1]) - 2.5*np.log10(x)) for x in fluxlc]
            absmag =[(np.min(point) - 5*(np.log10(dispc.value) - 1)) for point in obsmag]
            
            fluxdatapoints = np.array(y[4][:,4], dtype = float)
            fluxdatapointsall = np.array(hml[:,4], dtype = float)

            fluxyerrall = np.array(hml[:,5], dtype = float)
            fluxyerr = np.array(y[4][:,5], dtype = float)
            
            
            obsdatapoints =np.array(y[4][:,2], dtype = float)
            absdatapoints = [(np.min(point) - 5*(np.log10(dispc.value) - 1)) for point in obsdatapoints ]
            obsmagyerr = np.array(y[4][:,3], dtype = float)
            absmagvalue =np.min(absmag)
            
            
            '''
            plt.subplot(6,6,jo)
            plt.plot(timearray, fluxlc, color= 'blue')
            plt.errorbar(np.array(hml[:,1], dtype = float), fluxdatapointsall , yerr=fluxyerrall, fmt = 'o', color= 'green')
            plt.errorbar(np.array(y[4][:,1], dtype = float), fluxdatapoints , yerr=fluxyerr, fmt = 'o', color= 'blue')
            plt.xlabel('Time /days')
            plt.ylabel('Flux') #Absolute Magnitude')
            plt.title(sn + ' absmag: ' + str(absmagvalue) + ' redchisq: ' + str(y[2]))
           # plt.tight_layout()
           # plt.gca().invert_yaxis()
            '''
            jo+=1
            


            if y[2] < redchilimit:
                if sntype[1] =='b':
                    peakabsmagvalueb.append(absmagvalue)
                else:
                    peakabsmagvaluec.append(absmagvalue)

                print sn, 'abs peak mag', absmagvalue
            
                
            

print jo
#plt.show()


#print j

    
plt.subplot(211)  
peakabsmagvalue = [peakabsmagvalueb,peakabsmagvaluec] 
plt.hist(peakabsmagvalue, bins = redchilimit/0.8, label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'], stacked=True)
#plt.hist(peakabsmagvaluec, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green',stacked=True)
plt.title("Peak Absolute Magnitude Histogram of TypeIb and TypeIc using the nugent-sn1bc model")
plt.xlabel("Absolute Magnitude")
plt.ylabel("Frequency")
plt.legend()

plt.subplot(212)            
redchi = [redchib,redchic]
plt.hist(redchi, bins = redchilimit/1, label = [str(ib) + ' Ib datapoints',str(ic) + ' Ic datapoints' ], color = ['blue','green'], stacked=True)
#plt.hist(redchic, bins = redchilimit/0.5, label = str(ic) + ' datapoints', color = 'green')

plt.title("Reduced Chi Squared Distribution of TypeIb/c using the nugent-sn1bc model")
plt.xlabel("Reduced Chi Squared")
plt.ylabel("Frequency")
plt.legend()
plt.show()



'''
f = open('/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/snlistmjd.txt', 'w')
for item in datalistmjd:
    for i in item:
        print i
        f.write(str(i) +' ')
    f.write("\n")
    
f.close()
'''