import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
import urllib as url
import math as math
from usno import usno 
from cmath import cos, sin, phase
from math import radians, degrees

fits1 = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'
s1 = pf.open(fits1)
# Read position from the FITS file and convert RA/DEC to degrees
# be sure to check that the header data is reliable. If not
# edit the position by hand.
ras = s1[0].header['ra']
des = s1[0].header['dec']
radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
dsgn = np.sign(float(des[0:2]))
dedeg = float(des[0:2]) + dsgn*float(des[4:5])/60. + dsgn*float(des[7:])/3600.
fovam = 3.0 # size of square search field in arc min
epoch = s1[0].header['equinoxu']
name,rad,ded,rmag = usno(radeg,dedeg,fovam,epoch)
print(rad)

plt.figure(1)
w = np.where(rmag < 28)[0]
plt.subplot(211)
plt.plot(rad[w],ded[w],'g.')
plt.locator_params(axis='x',nbins=4)
plt.locator_params(axis='y',nbins=4)
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]')
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
#plt.xlim([106.0,105.0]) # reverse the x-axis direction


rar0 = radians((np.amax(rad[w])+np.amin(rad[w]))/2) #RA0 in radian
der0 = radians((np.amax(ded[w])+np.amin(ded[w]))/2) #Declination0 in radian

rar = []
der = []
X = []
Y = []
x=[]
x1=[]
y=[]
y1=[]
for i in w:
    rar.append(radians(rad[i])) #RA in radian
    der.append(radians(ded[i])) #Declination in radian

j = 0
while (j < len(rar)):
    X.append(-((cos(der[j])*sin(rar[j]-rar0))/(cos(der0)*cos(der[j])*cos(rar[j]-rar0)+sin(der[j])*sin(der0))))
    Y.append(-((sin(der0)*cos(der[j])*cos(rar[j]-rar0)-cos(der0)*sin(der[j]))/(cos(der0)*cos(der[j])*cos(rar[j]-rar0)+sin(der[j])*sin(der0))))
    j+=1


k = 0
while (k < len(X)):
    x.append(abs(30.84*(X[k]/0.000015)))
    y.append(abs(30.84*(Y[k]/0.000015)))        
    k+=1
  
l=0
while (l < len(x)):
    x1.append(x[l].real)
    y1.append(y[l].real)
    l+=1
    
print(x1)
print(y1)

plt.subplot(212)
plt.plot(x1, y1, 'r.')
plt.xlabel('x [Pixel]')
plt.ylabel('y [Pixel]')
plt.xlim((0,1000))
plt.ylim((0,1000))
plt.show()
