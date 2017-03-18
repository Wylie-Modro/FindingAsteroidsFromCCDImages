import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub
import urllib as url
import math as math
from usno import usno 
from cmath import cos, sin, phase
from math import radians, degrees
from FindStars import LocateMainPeakRanges
from compareDataUsno import GetCCDxyFromUSNOFits

print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

xRotated, yRotated, XCyl, YCyl = GetCCDxyFromUSNOFits(filename)

fp = 190020

i=0
matrixB = []
while i < len(XCyl):
    matrixB.append([fp*XCyl[i], fp*YCyl[i], 1])
    i+=1
    
print(matrixB)
 
c =   