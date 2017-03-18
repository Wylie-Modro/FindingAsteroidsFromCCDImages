import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub
import urllib as url
import math as math
import numpy.linalg 
from usno import usno 
from cmath import cos, sin, phase
from math import radians, degrees
from FindStars import LocateMainPeakRanges
from compareDataUsno import GetCCDxyFromUSNOFits
from numpy.linalg import inv, det

print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

xRotated, yRotated, XCyl, YCyl = GetCCDxyFromUSNOFits(filename)

fp = 190020

i=0
matrixB = []
while i < len(XCyl):
    matrixB.append([fp*XCyl[i], fp*YCyl[i], 1])
    i+=1

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

affineTransformation = np.matrix([[fp*2.85341734, -fp*0.76457087, 383.12313296],
                                 [fp*0.76457087, fp*2.85341734, 556.1085011], [0,0,1]])
#print((np.linalg.det(affineTransformation))**(1/2))
print("affineTransformation: " + str(affineTransformation))
print("affineTransformation.shape: " + str(affineTransformation.shape))

leastSquareX = np.matrix.dot(np.transpose((xRotatedMatrix - np.matrix.dot(matrixB, cX))),(xRotatedMatrix - np.matrix.dot(matrixB, cX)))
leastSquareY = np.matrix.dot(np.transpose((yRotatedMatrix - np.matrix.dot(matrixB, cY))),(yRotatedMatrix - np.matrix.dot(matrixB, cY)))
print(leastSquareX)
print(leastSquareY)


filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

x = pf.getdata(filename)
hdr = pf.getheader(filename)
xb = bsub(x,hdr.get('cover'))


flat = pf.getdata(filename)
fhdr = pf.getheader(filename)
flatb = bsub(flat,hdr.get('cover')) # Bias subtract
flatb = flatb/np.median(flatb) # normalize


starRows, starCols = LocateMainPeakRanges(flatb, 10., 5.0, 25, 25) 

littlexMatrix = []
i = 0
while (i < len(starRows)):
    littlexMatrix.append([starRows[i], starCols[i], 1.])
    i +=1
littlexMatrix = np.matrix(littlexMatrix)
print("littlexMatrix: " + str(littlexMatrix))  
bigXMatrix  = np.dot(inv(affineTransformation), np.transpose(littlexMatrix))
print("bigXMatrix: " + str(bigXMatrix))
print("bigXMatrix.shape: " + str(bigXMatrix.shape))

bigXValues = bigXMatrix[0]
bigYValues = bigXMatrix[1]
bigOneValues = bigXMatrix[2]



result1 = list(np.array(bigXValues).reshape(-1,))
bigXValuesList = list(np.array(result1).reshape(-1,))

result2 = list(np.array(bigYValues).reshape(-1,))
bigYValuesList = list(np.array(result2).reshape(-1,))


starNMinDistance = []
tempCombinedDiff = []
p = 0

print("length",len(XCyl))
while (p < len(bigXValuesList)):
    tempCombinedDiff[:] = [] #clears list
    w=0
    while (w < len(XCyl)): 
        print(bigXValuesList[p])
        print(XCyl[w])
        print(math.sqrt(((bigXValuesList[p] - XCyl[w])**2)+((bigYValuesList[p] - YCyl[w])**2)))
        tempCombinedDiff.append(math.sqrt(((bigXValuesList[p] - XCyl[w])**2)+((bigYValuesList[p] - YCyl[w])**2)))
        print(tempCombinedDiff)
        print(p,w)
        w+=1
    starNMinDistance.append([bigXValuesList[p], bigYValuesList[p], min(tempCombinedDiff)])
    p+=1


print("starXNMinDistance: " + str(starNMinDistance))


XPlottingOnXAxis = [i[0] for i in starNMinDistance]
XplottingOnYAxis = [i[2] for i in starNMinDistance]

yPlottingOnXAxis = [i[1] for i in starNMinDistance]
yPlottingOnYAxis = [i[2] for i in starNMinDistance]

X_Differences, = plt.plot(XPlottingOnXAxis, XplottingOnYAxis, "g^", label = 'x')
Y_Differences, = plt.plot(yPlottingOnXAxis, yPlottingOnYAxis, "r^", label = 'y')

plt.legend(handles=[X_Differences, Y_Differences])
plt.title("Residual Errors in the Standard Coordinates")
plt.xlabel("X or Y [standard coordinates]")
plt.ylabel("Error [standard coordinates]")
#plt.xlim((0,1000))
#plt.ylim((0,200))
plt.show()





