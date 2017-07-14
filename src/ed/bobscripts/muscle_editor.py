'''
Created on 28 Jun 2016

@author: mate

Take in a muscle alignment from http://www.ebi.ac.uk/Tools/msa/ and concatenate the lines properly to form one single line
'''

with open("/home/mate/workspace/katamari/src/root/ed/datafiles/muscle-I20160628-172926-0761-43925401-pg.clw","r") as inpF:

  sL = []
  
  lineCount = 0
  startFlag = True
  for inpLine in inpF:
    inpLine = inpLine[33:].strip(" \n")
    print repr(inpLine)
    if startFlag:
      startFlag = False
      continue
    inpS = inpLine.rstrip()
    if inpS == "": 
      lineCount = 0
      print "+++"
      continue
    try:
      sL[lineCount] += inpS 
      print len(inpS)
    except IndexError:
      sL.append(inpS)
    lineCount += 1
    if lineCount == 7: break

with open("/home/mate/workspace/katamari/src/root/ed/datafiles/muscle-I20160628-172926-0761-43925401-pg.txt","w") as outF:
  for sI in sL:
    outF.write(sI+"\n")
    