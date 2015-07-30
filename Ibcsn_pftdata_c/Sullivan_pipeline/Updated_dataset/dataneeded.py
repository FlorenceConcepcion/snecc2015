import numpy as np
import subprocess, math, sncosmo, astropy, os
import matplotlib.pyplot as plt


f = np.loadtxt('/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/Sullivan_pipeline/completesnlist.txt', dtype= str,skiprows=2,  usecols = (0,5,6), delimiter = '|')

z = []
liste = {}
for i in f:
    liste[(i[0].strip(' '))] = float(i[2])
    #print i[0].strip(' '), i[1].strip('  '), float(i[2])
print 'liste', liste
print len(liste)
completesnlist = ['09awk', '09dfk', '09dh', '09dha', '09dzt', '09fsr', '09iqd', '09ps', '09q', '09sk', '09ut', '10aavz', '10abck', '10acbu', '10acff', '10acgq', '10bhu', '10bip', '10bzf', '10ciw', '10cs', '10eqi', '10fbv', '10feq', '10fia', '10fmx', '10gvb', '10hfe', '10hie', '10inj', '10iue', '10kui', '10lbo', '10ood', '10osn', '10pbi', '10qif', '10qqd', '10qts', '10svt', '10tqi', '10tqv', '10vgv', '10vnv', '10wal', '10xem', '10xik', '10xjr', '10yow', '10ysd', '10zcn', '11awe', '11bli', '11bov', '11cmh', '11gcj', '11hyg', '11ilr', '11img', '11ixk', '11izq', '11jgj', '11kaa', '11klg', '11lbm', '11lmn', '11mnb', '11mwk', '11pnq', '11qcj', '11qiq', '11rfh', '11rka', '12as', '12bwq', '12cde', '12cjy', '12dcp', '12dtf', '12eaw', '12elh', '12fes', '12fgw', '12gdy', '12gps', '12grr', '12gty', '12gvr', '12gzk', '12hni', '12hvv', '12jxd', '12ktu', '12ldy', '12lpo', '12ltw', '12lvt', '12mfx']


print len(completesnlist)

#creating list of available supernova data
datadir = '/home/fcm1g13/Documents/Supernova/Ibcsn_pftdata_c/Sullivan_pipeline/data_files'

cclist = []; snlist = []
for i in os.walk(datadir):
    for j in i[2]:
        if j[-1:] == 'R':
            j[3:j.index('.')]
            cclist.append(j[3:j.index('.')])
#print cclist
print len(cclist)



snsofar = []
        
for i in completesnlist:
    try:
        cclist.index(i)
        snsofar.append(i)
    except:
        print i        
        
print snsofar