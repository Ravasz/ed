'''
Created on 21 Sep 2017

@author: mate

Open .csv files made from western blot quantifications by the image studio software. 
These come from the report function in the western analysis section of the full IS suite. 
The xls files exported there where then converted in bulk into .csv using Automator on Mac.

Extract OST bands, OT1 bands and CSK bands. Count number of bands and write to file a reformatted .csv dataset for graphing in R.
'''

import os, numpy

cskSizeList = [] # stores CSK band sizes in KDa
cskSignalList = [] # CSK signal intensities

ostSizeList = [[]] # list of list PTPN22 band sizes. First list item is the largest protein in the lane, second is the second largest in the lane and so on
ostSignalList = [[]] # list of PTPN22 signal intensities. not used for anything for the time being

ot1SizeList = [[]]
ot1SignalList = [[]]




directoryS = os.path.join(os.path.split(os.path.dirname(__file__))[0], "data", "size_measurements")

outF = open(os.path.join(directoryS , "ost-ptpn22-sizes-full.csv"),"w")

for fileName in os.listdir(directoryS): # go through all files in the folder
    if fileName.endswith(".csv") and "2017-09" not in fileName: 
        # print(os.path.join(directoryS, fileName))
        inpFile = open(os.path.join(directoryS, fileName),"rU")
        headerFlag = True
        
        for inpLine in inpFile:
          if headerFlag:
            headerFlag = False
            continue
          # print inpLine
          inpList = inpLine.rstrip().split(",")
          
          if len(inpList) < 6: continue # get rid of bad lines
          
          if inpList[1] == "800" and "csk" in inpList[3]: # store Csk signal here
            #print inpList
            cskSizeList.append(float(inpList[5]))
            cskSignalList.append(float(inpList[4]))
          
          if inpList[1] == "700" and "ost" in inpList[3]: # store ost sizes and signal intensities
            if "-" in inpList[3]:
              ostPos = int(inpList[3].split("-")[-1]) # bands are labelled ost-1, ost-2 etc. this records which band this is
              try:
                ostSizeList[ostPos-1].append(float(inpList[5]))
              except IndexError:
                ostSizeList.append([])
                ostSizeList[ostPos-1].append(float(inpList[5]))
              try:
                ostSignalList[ostPos-1].append(float(inpList[4]))
              except IndexError:
                ostSignalList.append([])
                ostSignalList[ostPos-1].append(float(inpList[4]))
                            
            else:
              ostSizeList[0].append(float(inpList[5]))
              ostSignalList[0].append(float(inpList[4]))

          if inpList[1] == "700" and "ot" in inpList[3]: # store ot-1 sizes and signal intensities
            if "-" in inpList[3]:
              ostPos = int(inpList[3].split("-")[-1]) # bands are labelled ost-1, ost-2 etc. this records which band this is
              try:
                ostSizeList[ostPos-1].append(float(inpList[5]))
              except IndexError:
                ostSizeList.append([])
                ostSizeList[ostPos-1].append(float(inpList[5]))
              try:
                ostSignalList[ostPos-1].append(float(inpList[4]))
              except IndexError:
                ostSignalList.append([])
                ostSignalList[ostPos-1].append(float(inpList[4]))
                            
            else:
              ostSizeList[0].append(float(inpList[5]))
              ostSignalList[0].append(float(inpList[4]))            
        
        continue
    else:
        continue

outDataF = open(os.path.join(directoryS , "ost-ptpn22-sizes-summary.csv"),"w")
listCount = 0
ostSingleList = []
print "ost bands:"
print ostSizeList
for ostSubList in ostSizeList:
  listCount += 1
  outFSegmented = open(os.path.join(directoryS , "ost-ptpn22-sizes" + str(listCount) + ".csv"),"w")
  for ostSubItem in ostSubList:
    ostSingleList.append(ostSubItem)
    outFSegmented.write(str(ostSubItem) + "\n")
  outFSegmented.close()
  
  print listCount

  print len(ostSubList)
  print round(sum(ostSubList) / float(len(ostSubList)),1)
  print round(numpy.std(ostSubList),1)
  print "---"

ostSingleList.sort()
for ostSize in ostSingleList:
  outF.write(str(ostSize) + "\n")

outF.close()

print "Csk: bands found, average size, SD:"
print len(cskSizeList)
print sum(cskSizeList) / float(len(cskSizeList))
print numpy.std(cskSizeList)
