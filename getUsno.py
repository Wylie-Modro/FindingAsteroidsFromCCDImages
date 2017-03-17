import astropy.io.fits as pf
import matplotlib.pyplot as plt
import numpy as np
import urllib as url
from usno import usno 

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
w = np.where(rmag < 17)[0]
plt.plot(rad[w],ded[w],'g.')
plt.locator_params(axis='x',nbins=4)
plt.locator_params(axis='y',nbins=4)
plt.tick_params('x',pad=10)
plt.ylabel('Dec [Deg]')
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
'plt.xlim([266.11,266.03]) # reverse the x-axis direction'
plt.show()
