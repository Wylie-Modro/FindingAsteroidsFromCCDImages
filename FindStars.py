import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub

def LocatePeakRanges(flatb, gap, threashhold):
        peakRanges=[]
        intensitiesDict = {}
        rowNum = 0
        for entry in flatb:
            colNum = 0
            rowNum +=1
            for subEntry in entry:
                colNum +=1
                #print("rowNum: " + str(rowNum), " colNum: " + str(colNum))
                colRowHack = int(str(rowNum)+str(colNum))
                print("colRowHack : " + str(colRowHack))
                #subEntryNpArray = np.ndarray(subEntry)
                #print(subEntry)
                intensitiesDict[colRowHack] = subEntry
        #print(intensitiesDict)  
        for key in intensitiesDict.keys():
            print("key: " + str(key))
            if key+gap >= 10241024:
                gap = float(int(gap/2))
            startValue = intensitiesDict[key]
            midValue = intensitiesDict[key+int(gap/2)]
            endValue = intensitiesDict[key+gap]
            if midValue - startValue > threashhold and midValue - endValue > threashhold :
                #print('Peak at:' + str(key + gap/2) +'!!!')
                peakRanges.extend([key, key + gap/2, key + gap])
            else:
                #print('No peak at:' + str(key + gap/2)+' intensitiesDict[key+gap]: ' + str(intensitiesDict[key+gap]))
                pass
        return peakRanges
        

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

peakRanges = LocatePeakRanges(flatb, 10., 0.2)
print(peakRanges)
#print(IsolatePeaks(peakRanges))


