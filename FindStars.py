import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from bsub import bsub


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
                #print("rowNum: " + str(rowNum), " colNum: " + str(colNum)
                #print("colRowHack : " + str(colRowHack))
                #subEntryNpArray = np.ndarray(subEntry)
                #print(subEntry)
                intensitiesDict[num] = [rowNum, colNum, subEntry]
                #print(intensitiesDict[num])
        #print(intensitiesDict)  
        
        draftLocationOfStarRanges = LocateRowPeakRanges(intensitiesDict, intensityGap, intensityThreshhold)
        print("draftLocationOfStarRanges: " + str(draftLocationOfStarRanges))
        print("len(draftLocationOfStarRanges): " + str(len(draftLocationOfStarRanges)))
        #print("intensitiesDict[test]"  + str(intensitiesDict[draftLocationOfPeakRanges[0]]) + str(intensitiesDict[draftLocationOfPeakRanges[1]]) + str(intensitiesDict[draftLocationOfPeakRanges[2]]))
        #confirmedLocationOfPeakRanges = LocateColPeakRanges(draftLocationOfPeakRanges, intensitiesDict, colThreshold, rowThreshold)
        groupedLocationOfStarRanges = IsolateStars(draftLocationOfStarRanges, intensitiesDict, colThreshold, rowThreshold)
        print("groupedLocationOfStasrRanges: " + str(groupedLocationOfStarRanges))
        print("len(groupedLocationOfStarRanges): " + str(len(groupedLocationOfStarRanges)))
        averagedStarLocations = GetAveStarLocFromAverage(groupedLocationOfStarRanges, intensitiesDict)
        print("averagedStarLocations: " + str(averagedStarLocations))
        print("len(averagedStarLocations): " + str(len(averagedStarLocations)))
        
        starRows = []
        starCols = []
        #starsRowsCols = np.ndarray(shape=(len(confirmedLocationOfPeakRanges), len(confirmedLocationOfPeakRanges)), dtype=int)
        for star in averagedStarLocations:
            #print(star)
            starRows.append(star[0])
            starCols.append(star[1])
            #np.append(starsRowsCols,[(intensitiesDict[i])[0], (intensitiesDict[i])[1]])
        print("starRows: " + str(starRows))
        print("starCols: " + str(starCols))
        plt.plot(starCols, starRows, "r.")
        plt.xlim([0.0, 1024.0])
        plt.ylim([0.0, 1024.0])
        plt.title("intensityGap: " + str(intensityGap) + "| intensityT " + str(intensityThreshhold) + "| colT: " + str(colThreshold) + "| rowT: " + str(rowThreshold))
        #plt.xlim([460.0, 464.0])
        #plt.ylim([435.0, 455.0])
        plt.show()
        
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
                #locationRowOfPeakRanges.extend([int(key), int(key + gap/2), int(key + gap)])
                locationRowOfPeakRanges.append(int(key + gap/2))
            else:
                #print('No peak at:' + str(key + gap/2)+' intensitiesDict[key+gap]: ' + str(intensitiesDict[key+gap]))
                pass
    return locationRowOfPeakRanges
        
'''        
def LocateColPeakRanges(locationRowOfPeakRanges, intensitiesDict, colThreshold, rowThreshold):
    # get col and row of each
    #for keyOuterLoop, valueOuterLoop in intensitiesDict.iteritems():
    potentialStarGroup = set()
    for peakNumOuterLoop in locationRowOfPeakRanges: 
        
        #print("peakNumOuterLoop: " + str(peakNumOuterLoop))
        rowOuterLoop = (intensitiesDict[peakNumOuterLoop])[0]
        #print("rowOuterLoop: " +str(rowOuterLoop))
        colOuterLoop = (intensitiesDict[peakNumOuterLoop])[1]
        intensityOuterLoop = (intensitiesDict[peakNumOuterLoop])[2]
        
        #for keyInnerLoop, valueInnerLoop in intensitiesDict.iteritems():
        for peakNumInnerLoop in locationRowOfPeakRanges:
            
            #print("peakNumInnerLoop: " + str(peakNumInnerLoop))
            rowInnerLoop = (intensitiesDict[peakNumInnerLoop])[0]
            #print("rowInnerLoop: " +str(rowInnerLoop))
            colInnerLoop = (intensitiesDict[peakNumInnerLoop])[1]
            intensityInnerLoop = (intensitiesDict[peakNumInnerLoop])[2]
            
            if (rowOuterLoop != rowInnerLoop and abs(rowOuterLoop - rowInnerLoop) < rowThreshold):
                if abs(colOuterLoop - colInnerLoop) < colThreshold:
                    if (colOuterLoop <= 275 and colOuterLoop >= 245): #Remove Line
                        pass
                    else: 
                        #potentialStars.add([peakNumOuterLoop, peakNumInnerLoop]) 
                        potentialStarGroup.add(peakNumOuterLoop)
                        potentialStarGroup.add(peakNumInnerLoop)
                else:
                    #add to new set
                    pass
    
    return potentialStarGroup
    
    # compare col with another, if close and if row is close and add range to another list  
    
''' 

def IsolateStars(locationRowOfPeakRanges, intensitiesDict, colThreshold, rowThreshold):
        peakSet = set()
        masterList = []
        tempList = []
        
        for peakNumOuterLoop in locationRowOfPeakRanges:
            rowO = (intensitiesDict[peakNumOuterLoop])[0]
            colO = (intensitiesDict[peakNumOuterLoop])[1]
            print("peakNumOuterLoop: " + str(peakNumOuterLoop))
            
            if (colO <= 275 and colO >= 245):
                locationRowOfPeakRanges.remove(peakNumOuterLoop)
                print("Deleted peakNumOuterLoop as was in line")
            else:  
            
                #print("----------------------------------")
                for peakNumInnerLoop in locationRowOfPeakRanges:
                    rowI = (intensitiesDict[peakNumInnerLoop])[0]
                    colI = (intensitiesDict[peakNumInnerLoop])[1]
                    print("peakNumInnerLoop: " + str(peakNumInnerLoop))
                    
                    #nextRow = (intensitiesDict[locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]])[0]
                    print("rowO: " + str(rowO))
                    print("rowI: " + str(rowI))
                    print("colO: " + str(colO))
                    print("colI: " + str(colI))
                    #print("nextRow: " + str(nextRow))
                    
                    #nextCol = (intensitiesDict[locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]])[1]
    
                    #print("nextCol: " + str(nextCol))
                    if (colI <= 275 and colI >= 245):
                            locationRowOfPeakRanges.remove(peakNumInnerLoop)
                            print("Deleted peakNumInnerLoop as was in line")
                    elif np.abs(rowO - rowI) <= rowThreshold and np.abs(colO - colI) <= colThreshold:
                        print("Through Threshhold")
                        if (rowO == rowI and colO == colI):
                            pass
                        else:
                            #print("Added peakNumOuterLoop: " + str(peakNumOuterLoop))
                            peakSet.add(peakNumInnerLoop)
                            
                            locationRowOfPeakRanges.remove(peakNumInnerLoop)
                            print("add: len(locationRowOfPeakRanges):" + str(len(locationRowOfPeakRanges)))
                            
                            #print("Added locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]: " + str(locationRowOfPeakRanges[locationRowOfPeakRanges.index(peakNumOuterLoop)+1]))
                            #print("peakSet after add: " + str(peakSet))
                    else: 
                        print("Rejected")
                        pass
                    
                    # when itterated peakNumOuterLoop for all values of peakNumInnerLoop
                peakSet.add(peakNumOuterLoop)
                locationRowOfPeakRanges.remove(peakNumOuterLoop)
                print("end loop: len(locationRowOfPeakRanges):" + str(len(locationRowOfPeakRanges)))
                
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
            

def MergeGroupsThatAreClose(intensitiesDict, masterList, colThreshold, rowThreshold):
    finalGroups = []
    for group1 in masterList:
        for pixel1 in group1:
            
            row1 = (intensitiesDict[pixel1])[0]
            col1 = (intensitiesDict[pixel1])[1]
            
            for group2 in masterList:
                for pixel2 in group2:
                    
                    row2 = (intensitiesDict[pixel2])[0]
                    col2 = (intensitiesDict[pixel2])[1]
                    
                    if (abs(row1 - row2) < rowThreshold and abs(col1 - col2) < colThreshold) and group1 != group2:
                        pass
                        
                


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


peakRanges = LocateMainPeakRanges(flatb, 10., 10.0, 25, 25) 
#print("peakRanges: " + str(peakRanges))
#print(IsolatePeaks(peakRanges))


