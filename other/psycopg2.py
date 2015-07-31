
import subprocess, math, sncosmo, astropy, os
import numpy as np
import matplotlib.pyplot as plt

import psycopg2
conn = psycopg2.connect("host=srv01050.soton.ac.uk dbname=frohmaier user=frohmaier password=rates")
cur=conn.cursor()
    
datadir = '/home/fcm1g13/Documents/Supernova/CorecollapseLC/ccdata/'

cclist = []
for i in os.walk(datadir):
    for j in i:
        for k in j:
            if k[-4:] == '.npy':
                cclist.append(k)
def hi(sn):
    redshift={}

    cur.execute("SELECT * from specinfo where ptfname=%s;", (str(sn[:-7]),))
    m=cur.fetchone()
    cur.close()
    m=np.array(m)
    redshift[sn[:-7]] = m[5]

    print m[5]
    
hi('10gtn_LC.npy')