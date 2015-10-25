import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt
import psycopg2
from astropy.cosmology import FlatLambdaCDM
from astropy.time import Time

#creating list of available supernova data
datadir = '../data_files'

cclist = []; snlist = []
for i in os.walk(datadir):
    for j in i:
        for k in j:
            if k[-1:] == 'R':
                cclist.append(k)
            

testlist = cclist[:]
#print testlist

types = {'10osn': 'Ic', '10ciw': 'Ic-BL', '12gty': 'Ic', '10wal': 'Ic', '12mfx': 'Ib', '10vgv': 'Ic-BL', '10acff': 'Ib', '10xem': 'Ic-BL', '10aavz': 'Ic-BL', '11pnq': 'Ib/c', '12ktu': 'Ic', '10qqd': 'Ic', '10qts': 'Ic-BL', '11klg': 'Ic', '12gdy': 'Ic-BL', '11cmh': 'Ic-BL', '12cde': 'Ib/c', '10hfe': 'Ic', '11ixk': 'Ic', '11bov': 'Ic', '10vwg': 'Ic', '10ysd': 'Ic-BL', '10cs': 'Ic-BL', '10bip': 'Ic', '10fia': 'Ib', '11oni': 'Ib/c', '10kui': 'Ib', '12lpo': 'Ib/c', '10pbi': 'Ib', '11all': 'Ic', '10hcw': 'Ib', '12gvr': 'Ib/c', '10fbv': 'Ib', '11img': 'Ic-BL', '10zcn': 'Ic', '12ldy': 'Ibn', '11izq': 'Ib', '11lmn': 'Ib/c', '11kaa': 'Ib', '10acgq': 'Ib', '12fsu': 'Ic', '11gcj': 'Ib/c', '12dcp': 'Ibn', '10svt': 'Ic', '10vnv': 'Ib', '12eaw': 'Ib', '12as': 'Ic-BL', '12hni': 'Ic', '12gps': 'Ib', '10bhu': 'Ic', '12fgw': 'Ic', '11bli': 'Ib/c', '12elh': 'Ic', '11hyg': 'Ic', '10hie': 'Ic', '10xik': 'Ic', '10inj': 'Ib', '10lbo': 'Ic', '10yow': 'Ic', '10qif': 'Ib', '12jxd': 'Ic', '10acbu': 'Ic', '10uhf': 'Ic', '11lbm': 'Ic-BL', '10uhn': 'Ic', '11rka': 'Ic', '10abck': 'Ib', '10feq': 'Ib', '10bzf': 'Ic-BL', '10gvb': 'Ic-BL', '11ilr': 'Ib', '12ltw': 'Ib', '10xjr': 'Ib', '10fmx': 'Ic', '10iue': 'Ic-BL', '10hgi': 'Ib/c', '12dtf': 'Ic', '12fes': 'Ib', '12gzk': 'Ic', '12grr': 'Ic-BL', '10xfh': 'Ic', '11qcj': 'Ic', '11rfh': 'Ic', '11jgj': 'Ic', '10bfz': 'Ic-BL', '10tqv': 'Ic', '12cjy': 'Ic', '11awe': 'Ib', '11mnb': 'Ic', '10eqi': 'Ib', '10ood': 'Ic', '10tqi': 'Ic', '12lvt': 'Ib', '12bwq': 'Ib', '11dhf': 'Ib', '11mwk': 'Ib/c', '11qiq': 'Ib', '12hvv': 'Ic'}

z = []
for i in testlist:
    try:
        types[i[3:i.index('.')]]
        z.append(i[3:i.index('.')])
    except KeyError:
        pass
        #print i[3:i.index('.')]
        
#print z
nodata = []
correctsn = []
completesnlist = []
for i in types:
    completesnlist.append( i)
    try:
        z.index(i)
        correctsn.append(i)
    except ValueError:
        nodata.append( i)
        
for i in completesnlist:
    print types[i]

print correctsn

completesnlist = ['10osn', '10ciw', '12gty', '10wal', '12mfx', '10vgv', '10acff', '10xem', '10bhu', '11pnq', '12ktu', '12fgw', '10qts', '11klg', '12gdy', '11cmh', '12cde', '10hfe', '11ixk', '11bov', '10vwg', '10ysd', '10cs', '10bip', '10fia', '11oni', '10acgq', '12lpo', '12elh', '12cjy', '10hcw', '12gvr', '10fbv', '11img', '10zcn', '12ldy', '11izq', '11lmn', '11kaa', '10kui', '12fsu', '11gcj', '12dcp', '10svt', '10vnv', '12eaw', '12as', '12hni', '12gps', '10aavz', '10qqd', '11bli', '10pbi', '11hyg', '10hie', '10xik', '10inj', '10lbo', '10yow', '10qif', '12jxd', '10acbu', '10uhf', '11lbm', '10uhn', '11rka', '10abck', '10feq', '10bzf', '10gvb', '11ilr', '12ltw', '10xjr', '10fmx', '10iue', '10hgi', '12dtf', '12fes', '12gzk', '12grr', '10xfh', '11qcj', '11rfh', '11jgj', '10bfz', '10tqv', '11all', '11awe', '11mnb', '10eqi', '10ood', '10tqi', '12lvt', '12bwq', '11dhf', '11mwk', '11qiq', '12hvv']



nodata = ['10vgv', '10acff', '12gdy', '10vwg', '10cs', '10acgq', '12ldy', '11lmn', '11kaa', '12as', '12gps', '10aavz', '10pbi', '10yow', '10acbu', '11lbm', '10abck', '12ltw', '10xfh', '11awe', '10ood', '11dhf']

nodata = ['10acff', '10cs', '10acgq', '12ldy', '11lmn',  '12as', '10aavz', '10acbu', '11lbm', '10abck', '11awe']

