import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub

print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1002.fits'

x = pf.getdata(filename)
hdr = pf.getheader(filename)
xb = bsub(x,hdr.get('cover'))
print(xb)

flat = pf.getdata(filename)
fhdr = pf.getheader(filename)
flatb = bsub(flat,hdr.get('cover')) # Bias subtract
flatb = flatb/np.median(flatb) # normalize
plt.figure()
plt.subplot(121)
imp = plt.imshow(flatb,cmap='gray_r',vmin=0.8,vmax=1.2)
plt.colorbar()
plt.subplot(122)
imp = plt.imshow(xb/(flatb),cmap='gray_r',vmin=50,vmax=500)
plt.colorbar()
plt.show()