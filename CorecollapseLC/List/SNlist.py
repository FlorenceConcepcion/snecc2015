
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

data=np.loadtxt('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist.txt', dtype = str)
file = open('/home/fcm1g13/Documents/Supernova/CorecollapseLC/List/SNlist reduced.txt', 'w')


for i in data:
    print i
    coordinates = ''.join(j + ' ' for j in i[1:7])
    c = SkyCoord(coordinates, unit=(u.hourangle, u.deg))
    file.write(str(i[0]) + ' ' +str(c.ra.degree) + ' ' + str(c.dec.degree) + '\n')
file.close()
