import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub
import urllib as url
import math as math
from usno import usno 
from cmath import cos, sin, phase
from math import radians, degrees

def LocateMainPeakRanges(flatb, intensityGap, intensityThreshhold, colThreshold, rowThreshold):
        intensitiesDict = {}
        rowNum = 0
        num = 0
        for entry in flatb:
            colNum = 0
            rowNum +=1
            for subEntry in entry:
                colNum += 1
                num += 1
                intensitiesDict[num] = [rowNum, colNum, subEntry]
        
        draftLocationOfStarRanges = LocatePeakRanges(intensitiesDict, intensityGap, intensityThreshhold)
        print("draftLocationOfStarRanges: " + str(draftLocationOfStarRanges))
        print("len(draftLocationOfStarRanges): " + str(len(draftLocationOfStarRanges)))
        groupedLocationOfStarRanges = IsolateStarPointsIntoGroups(draftLocationOfStarRanges, intensitiesDict, colThreshold, rowThreshold)
        print("groupedLocationOfStasrRanges: " + str(groupedLocationOfStarRanges))
        print("len(groupedLocationOfStarRanges): " + str(len(groupedLocationOfStarRanges)))
        averagedStarLocations = GetAveStarLocFromAverage(groupedLocationOfStarRanges, intensitiesDict)
        print("averagedStarLocations: " + str(averagedStarLocations))
        print("len(averagedStarLocations): " + str(len(averagedStarLocations)))
        
        starRows = []
        starCols = []
        for star in averagedStarLocations:
            #print(star)
            starRows.append(star[0])
            starCols.append(star[1])
        print("starRows: " + str(starRows))
        print("starCols: " + str(starCols))
        plt.plot(starCols, starRows, "r.")
        plt.xlim([0.0, 1024.0])
        plt.ylim([0.0, 1024.0])
        plt.title("intensityGap: " + str(intensityGap) + "| intensityT " + str(intensityThreshhold) + "| colT: " + str(colThreshold) + "| rowT: " + str(rowThreshold))
        
def LocatePeakRanges(intensitiesDict, gap, threshhold):
    locationRowOfPeakRanges = []
    for key, value in intensitiesDict.iteritems():
            row = value[0]
            
            if row+gap >= 1024:
                gap = float(int(gap/2))
            startValue = (intensitiesDict[key])[2]
            midValue = (intensitiesDict[key+int(gap/2)])[2]
            endValue = (intensitiesDict[key+gap])[2]
            if midValue - startValue > threshhold and midValue - endValue > threshhold :
                locationRowOfPeakRanges.append(int(key + gap/2))
            else:
                pass
    return locationRowOfPeakRanges
        

def IsolateStarPointsIntoGroups(locationRowOfPeakRanges, intensitiesDict, colThreshold, rowThreshold):
        peakSet = set()
        masterList = []
        tempList = []
        
        for peakNumOuterLoop in locationRowOfPeakRanges:
            rowO = (intensitiesDict[peakNumOuterLoop])[0]
            colO = (intensitiesDict[peakNumOuterLoop])[1]
            #print("peakNumOuterLoop: " + str(peakNumOuterLoop))
            
            if (colO <= 275 and colO >= 245):
                locationRowOfPeakRanges.remove(peakNumOuterLoop)
                #print("Deleted peakNumOuterLoop as was in line")
            else:  
            
                for peakNumInnerLoop in locationRowOfPeakRanges:
                    rowI = (intensitiesDict[peakNumInnerLoop])[0]
                    colI = (intensitiesDict[peakNumInnerLoop])[1]
                    #print("peakNumInnerLoop: " + str(peakNumInnerLoop))
                    
                    #print("rowO: " + str(rowO))
                    #print("rowI: " + str(rowI))
                    #print("colO: " + str(colO))
                    #print("colI: " + str(colI))
                    
                    if (colI <= 275 and colI >= 245):
                            locationRowOfPeakRanges.remove(peakNumInnerLoop)
                            #print("Deleted peakNumInnerLoop as was in line")
                    elif np.abs(rowO - rowI) <= rowThreshold and np.abs(colO - colI) <= colThreshold:
                        #print("Through Threshhold")
                        if (rowO == rowI and colO == colI):
                            pass
                        else:
                            #print("Added peakNumOuterLoop: " + str(peakNumOuterLoop))
                            peakSet.add(peakNumInnerLoop)
                            
                            locationRowOfPeakRanges.remove(peakNumInnerLoop)
                            #print("add: len(locationRowOfPeakRanges):" + str(len(locationRowOfPeakRanges)))
                            
                            #print("Added locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]: " + str(locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]))
                            #print("peakSet after add: " + str(peakSet))
                    else: 
                        #print("Rejected")
                        pass
                    
                    # when itterated peakNumOuterLoop for all values of peakNumInnerLoop
                peakSet.add(peakNumOuterLoop)
                locationRowOfPeakRanges.remove(peakNumOuterLoop)
                #print("end loop: len(locationRowOfPeakRanges):" + str(len(locationRowOfPeakRanges)))
                
                if len(peakSet) > 0:
                    for peak in peakSet:
                        #print("peakSet: " + str(peakSet))
                        tempList.append(peak)
                    copyOftempList = list(tempList)
                    #print("tempList: "  + str(tempList))
                    masterList.append(copyOftempList)
                    #print("masterList: " + str(masterList))
                    del tempList[:]
                    peakSet.clear()
                
                #lastRow = (intensitiesDict[locationRowOfPeakRanges[-1]])[0]
                #lastCol = (intensitiesDict[locationRowOfPeakRanges[-1]])[1]
                #print("lastRow: " + str(lastRow))
                #print("lastCol: " + str(lastCol))
                '''
                if nextRow == lastRow and nextCol == lastCol:
                    for peak in peakSet:
                        tempList.append(peak)
                    masterList.append(tempList)
                '''
                
        print("About to exit")
        return masterList    
               # print("----------------------------------")
            


def GetAveStarLocFromAverage(masterList, intensitiesDict):
    listOfAveragedStarPositions = []
    for group in masterList:
        #print("group: " + str(group))
        sumCol = 0
        sumRow = 0
        for value in group:
            #print("value: " + str(value))
            sumRow += (intensitiesDict[value])[0]
            sumCol += (intensitiesDict[value])[1]
        if len(group) > 0:
            groupRowAverage = sumRow / len(group)
            groupColAverage = sumCol / len(group)
            listOfAveragedStarPositions.append([groupRowAverage, groupColAverage])
    return listOfAveragedStarPositions

print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

x = pf.getdata(filename)
hdr = pf.getheader(filename)
xb = bsub(x,hdr.get('cover'))


flat = pf.getdata(filename)
fhdr = pf.getheader(filename)
flatb = bsub(flat,hdr.get('cover')) # Bias subtract
flatb = flatb/np.median(flatb) # normalize


peakRanges = LocateMainPeakRanges(flatb, 10., 5.0, 25, 25) 
#print("peakRanges: " + str(peakRanges))
#print(IsolatePeaks(peakRanges))

#print("peakRanges: " + str(peakRanges))
#print(IsolatePeaks(peakRanges))

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
'''
plt.subplot(211)
plt.plot(rad[w],ded[w],'g.')
plt.locator_params(axis='x',nbins=4)
plt.locator_params(axis='y',nbins=4)
plt.tick_params('x',pad=10)
plt.xlabel('RA [Deg]')
plt.ylabel('Dec [Deg]')
plt.ticklabel_format(useOffset=False)
plt.axis('scaled')
plt.xlim([106.0,105.0]) # reverse the x-axis direction
'''

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
xRotated=[]
yRotated=[]
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
    x.append(abs(20*(X[k]/0.000015)+280))
    y.append(abs(20*(Y[k]/0.000015)+670))        
    k+=1
  
l=0
while (l < len(x)):
    x1.append(x[l].real)
    y1.append(y[l].real)
    l+=1

m=0
thetaDegrees=31
thetaRadians=radians(thetaDegrees)
print(type(x1))
print(type(cos(thetaRadians)).real)
while (m < len(x1)):
    xRotated.append((x1[m]*(cos(thetaRadians)).real - y1[m]*(sin(thetaRadians)).real))
    yRotated.append((x1[m]*(sin(thetaRadians)).real + y1[m]*(sin(thetaRadians)).real))
    m+=1
    
print(xRotated)
print(yRotated)
#plt.subplot(212)
plt.plot(xRotated, yRotated, 'b.')
plt.xlabel('x [Pixel]')
plt.ylabel('y [Pixel]')
plt.xlim((0,1000))
plt.ylim((0,1000))
plt.show()

