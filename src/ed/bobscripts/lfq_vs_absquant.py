'''
Created on 22 Jun 2016

@author: mate
prepare .csv file for Rstudio to plot correlation between LFQ fold changes 
between KO and WT PTPN22 and abs quant fold changes calculated with the histone ruler method.
'''

from math import log

print "this is lfq vs absquant"

with open("/home/mate/workspace/katamari/src/ed/bob/processed/abs_quant_24H_T4-14-24-27-05-20164.txt","r") as inpF:
  with open("/home/mate/workspace/katamari/src/ed/bob/processed/abs_quant_24H_T4-14-24-27-05-20165.txt","w") as outF:
    headerFlag = True
    for inpLine in inpF:
      if headerFlag:
        headerFlag = False
        outF.write(inpLine.rstrip() + ",FC-LFQ,FC-absquant\n")
        print inpLine.rstrip().split(",")[13:19]
        continue
      inpL = inpLine.rstrip().split(",")
      fcLFQ = (float(inpL[16]) +float(inpL[17]) +float(inpL[18]) + 1) / (float(inpL[13]) + float(inpL[14]) + float(inpL[15]) + 1) # wtLFQ / koLFQ
      fcLFQ = log(fcLFQ,2)
      if fcLFQ<-10:
        print inpL[13:19]
      fcAbsQuant = (float(inpL[27]) + 1) / (float(inpL[26]) + 1) # avg abs quant wt / avg abs quant KO
      fcAbsQuant = log(fcAbsQuant,2)
      for inpI in inpL:
        outF.write(inpI)
        outF.write(",")
      outF.write(str(fcLFQ))
      outF.write(",")
      outF.write(str(fcAbsQuant))
      outF.write("\n")
      
      