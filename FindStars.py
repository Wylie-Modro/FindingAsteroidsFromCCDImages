import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub


def LocateMainPeakRanges(flatb, gap, threshhold):
        intensitiesDict = {}
        rowNum = 0
        num = 0
        for entry in flatb:
            colNum = 0
            rowNum +=1
            for subEntry in entry:
                colNum += 1
                num += 1
                #print("rowNum: " + str(rowNum), " colNum: " + str(colNum))
                colRowHack = int(str(rowNum)+str(colNum))
                #print("colRowHack : " + str(colRowHack))
                #subEntryNpArray = np.ndarray(subEntry)
                #print(subEntry)
                intensitiesDict[num] = [rowNum, colNum, subEntry]
                #print(intensitiesDict[num])
        #print(intensitiesDict)  
        
        draftLocationOfPeakRanges = LocateRowPeakRanges(intensitiesDict, gap, threshhold)
        print("draftLocationOfPeakRanges: " + str(draftLocationOfPeakRanges))
        confirmedLocationOfPeakRanges = LocateColPeakRanges(draftLocationOfPeakRanges, intensitiesDict, 5, 5)
        print("confirmedLocationOfPeakRanges: " + str(confirmedLocationOfPeakRanges))

     
def LocateRowPeakRanges(intensitiesDict, gap, threshhold):
    locationRowOfPeakRanges = []
    for key, value in intensitiesDict.iteritems():
            #print("key: " + str(key))
            
            row = value[0]
            col = value[1]
            intensity = value[2]
            
            if row+gap >= 1024:
                gap = float(int(gap/2))
            startValue = (intensitiesDict[key])[2]
            midValue = (intensitiesDict[key+int(gap/2)])[2]
            endValue = (intensitiesDict[key+gap])[2]
            if midValue - startValue > threshhold and midValue - endValue > threshhold :
                #print('Peak at:' + str(key + gap/2) +'!!!')
                locationRowOfPeakRanges.extend([key, key + gap/2, key + gap])
            else:
                #print('No peak at:' + str(key + gap/2)+' intensitiesDict[key+gap]: ' + str(intensitiesDict[key+gap]))
                pass
    return locationRowOfPeakRanges
        
        
def LocateColPeakRanges(locationRowOfPeakRanges, intensitiesDict, colThreshold, rowThreshold):
    # get col and row of each
    #for keyOuterLoop, valueOuterLoop in intensitiesDict.iteritems():
    potentialStars = []
    for peakNumOuterLoop in locationRowOfPeakRanges: 
        
        rowOuterLoop = (intensitiesDict[peakNumOuterLoop])[0]
        print("rowOuterLoop: " +str(rowOuterLoop))
        colOuterLoop = (intensitiesDict[peakNumOuterLoop])[1]
        intensityOuterLoop = (intensitiesDict[peakNumOuterLoop])[2]  
          
        #for keyInnerLoop, valueInnerLoop in intensitiesDict.iteritems():
        for peakNumInnerLoop in locationRowOfPeakRanges:
            
            rowInnerLoop = (intensitiesDict[peakNumInnerLoop])[0]
            print("rowInnerLoop: " +str(rowInnerLoop))
            colInnerLoop = (intensitiesDict[peakNumInnerLoop])[1]
            intensityInnerLoop = (intensitiesDict[peakNumInnerLoop])[2]
            
            if (rowOuterLoop != rowInnerLoop and abs(rowOuterLoop - rowInnerLoop) < rowThreshold):
                if abs(colOuterLoop - colInnerLoop) < colThreshold:
                    potentialStars.append([peakNumOuterLoop, peakNumInnerLoop]) 
    
    return potentialStars
    
    # compare col with another, if close and if row is close and add range to another list  
    
    

def IsolatePeaks(peakRanges):
        peakSet = set()
        masterList = []
        tempList = []
        for peak1 in peakRanges:
            if np.abs(peak1 - peakRanges[peakRanges.index(peak1)+1]) < 10:
                peakSet.add(peak1)
                peakSet.add(peakRanges[peakRanges.index(peak1)+1])
            else: 
                if len(peakSet) > 0:
                    for peak in peakSet:
                        tempList.append(peak)
                    copyOftempList = list(tempList)
                    masterList.append(copyOftempList)
                    del tempList[:]
                    peakSet.clear()
            
            if peakRanges[peakRanges.index(peak1)+1] == peakRanges[-1]:
                for peak in peakSet:
                    tempList.append(peak)
                masterList.append(tempList)
                return masterList      
        




print('Asteroids do not concern me, Admiral. - Darth Vader')

filename = 'data-2017-03-02-nickel-Shelley.Wright/d1060.fits'

x = pf.getdata(filename)
hdr = pf.getheader(filename)
xb = bsub(x,hdr.get('cover'))


flat = pf.getdata(filename)
fhdr = pf.getheader(filename)
flatb = bsub(flat,hdr.get('cover')) # Bias subtract
flatb = flatb/np.median(flatb) # normalize


print(flatb.size)

i = 0
for entry in flatb:
    i+=1
print(i)
print(i**2)

peakRanges = LocateMainPeakRanges(flatb, 10., 0.5)
print("peakRanges: " + str(peakRanges))
#print(IsolatePeaks(peakRanges))


