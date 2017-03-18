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
from numpy.linalg import inv

print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

xRotated, yRotated, XCyl, YCyl = GetCCDxyFromUSNOFits(filename)

fp = 190020

i=0
matrixB = []
while i < len(XCyl):
    matrixB.append([fp*XCyl[i], fp*YCyl[i], 1])
    i+=1

print("Under Pressure")
matrixA = [[1,2,3],[4,5,6],[7,8,9]]
matrixA = np.matrix(matrixA)
matrixAT = np.transpose(matrixA)
da = np.dot(matrixAT, matrixA)

xRotatedMatrix = np.transpose(np.matrix(xRotated))
yRotatedMatrix = np.transpose(np.matrix(yRotated))

matrixB = np.matrix(matrixB)

matrixBT = np.transpose(matrixB)
dpmatrixBTmatrixB = np.matrix(np.matrix.dot(np.matrix(matrixBT), np.matrix(matrixB)))
d = inv(dpmatrixBTmatrixB)
matrixBTdotxRotated = np.dot(np.transpose(matrixB),xRotatedMatrix)
matrixBTdotyRotated = np.dot(np.transpose(matrixB),yRotatedMatrix)
cX = np.dot(d, matrixBTdotxRotated)
print(cX)
cY = np.dot(d, matrixBTdotyRotated)
print(cY)



